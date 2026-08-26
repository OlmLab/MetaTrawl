"""DuckDB project database for MetaTrawl."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import os
import re
import time
from typing import Callable, Iterable

import duckdb
import polars as pl

from metatrawl import allele_mask


ExportProgressCallback = Callable[[dict[str, object]], None]
ACCESSION_PATTERN = re.compile(r"(GC[AF]_\d+(?:\.\d+)?)", re.IGNORECASE)


SCHEMA_SQL = """
CREATE TABLE IF NOT EXISTS sra_runs (
    run_id VARCHAR PRIMARY KEY,
    status VARCHAR NOT NULL DEFAULT 'added',
    added_at DOUBLE NOT NULL,
    updated_at DOUBLE NOT NULL,
    deleted_at DOUBLE
);

CREATE TABLE IF NOT EXISTS samples (
    sample_id VARCHAR PRIMARY KEY,
    run_id VARCHAR NOT NULL,
    status VARCHAR NOT NULL,
    created_at DOUBLE NOT NULL,
    updated_at DOUBLE NOT NULL
);

CREATE TABLE IF NOT EXISTS profiles (
    run_id VARCHAR PRIMARY KEY,
    profile_file VARCHAR NOT NULL,
    genome_stats_file VARCHAR NOT NULL,
    gene_stats_file VARCHAR,
    sylph_abundance_file VARCHAR NOT NULL,
    profile_storage_mode VARCHAR NOT NULL DEFAULT 'full',
    profile_min_cov INTEGER,
    created_at DOUBLE NOT NULL,
    updated_at DOUBLE NOT NULL
);

CREATE TABLE IF NOT EXISTS profile_positions (
    sample_id VARCHAR NOT NULL,
    chrom VARCHAR NOT NULL,
    pos BIGINT NOT NULL,
    genome VARCHAR NOT NULL,
    A USMALLINT NOT NULL,
    C USMALLINT NOT NULL,
    G USMALLINT NOT NULL,
    T USMALLINT NOT NULL,
    ref_base_bitmask UTINYINT
);

CREATE TABLE IF NOT EXISTS profile_storage (
    id UTINYINT PRIMARY KEY CHECK (id = 1),
    mode VARCHAR NOT NULL,
    format_version INTEGER NOT NULL,
    min_cov INTEGER,
    codec VARCHAR,
    allele_bit_order VARCHAR
);

CREATE TABLE IF NOT EXISTS allele_mask_reference_segments (
    segment_id UBIGINT PRIMARY KEY,
    genome VARCHAR NOT NULL,
    chrom VARCHAR NOT NULL,
    segment_ordinal INTEGER NOT NULL,
    start_pos BIGINT NOT NULL,
    span BIGINT NOT NULL,
    reference_mask BLOB NOT NULL,
    reference_hash VARCHAR NOT NULL,
    UNIQUE (genome, chrom)
);

CREATE TABLE IF NOT EXISTS allele_mask_profile_blocks (
    sample_id VARCHAR NOT NULL,
    segment_id UBIGINT NOT NULL,
    presence BLOB NOT NULL,
    deviation BLOB NOT NULL,
    covered_positions BIGINT NOT NULL,
    payload_hash VARCHAR NOT NULL,
    PRIMARY KEY (sample_id, segment_id)
);

CREATE TABLE IF NOT EXISTS allele_mask_migration_state (
    sample_id VARCHAR PRIMARY KEY,
    status VARCHAR NOT NULL,
    source_rows BIGINT,
    compressed_blocks BIGINT,
    covered_positions BIGINT,
    error VARCHAR,
    updated_at DOUBLE NOT NULL
);

CREATE TABLE IF NOT EXISTS genome_stats (
    sample_id VARCHAR NOT NULL,
    genome VARCHAR NOT NULL,
    coverage DOUBLE,
    breadth DOUBLE,
    ber DOUBLE,
    ref_ani DOUBLE,
    coverage_median DOUBLE,
    coverage_std DOUBLE,
    genome_length BIGINT,
    gap_mean DOUBLE,
    gap_std DOUBLE,
    "5x_cov_sites" BIGINT,
    heterogeneity DOUBLE,
    fug DOUBLE,
    reads_mapped BIGINT,
    conANI_reference DOUBLE,
    SNS_count BIGINT,
    SNV_count BIGINT,
    presence VARCHAR,
    genome_taxonomy VARCHAR
);

CREATE TABLE IF NOT EXISTS gene_stats (
    sample_id VARCHAR NOT NULL,
    genome VARCHAR,
    gene VARCHAR NOT NULL,
    coverage DOUBLE,
    breadth DOUBLE,
    ber DOUBLE,
    ref_ani DOUBLE,
    length BIGINT
);

CREATE TABLE IF NOT EXISTS sylph_abundance (
    sample_id VARCHAR NOT NULL,
    genome VARCHAR,
    accession VARCHAR,
    abundance DOUBLE
);

CREATE TABLE IF NOT EXISTS cache_genomes (
    accession VARCHAR PRIMARY KEY,
    genome_fasta VARCHAR NOT NULL,
    gene_fasta VARCHAR NOT NULL,
    status VARCHAR NOT NULL,
    updated_at DOUBLE NOT NULL
);

CREATE TABLE IF NOT EXISTS matrix_stores (
    matrix_id VARCHAR PRIMARY KEY,
    genome VARCHAR NOT NULL,
    matrix_file VARCHAR NOT NULL,
    profile_count BIGINT NOT NULL,
    storage_layout VARCHAR NOT NULL,
    min_coverage DOUBLE,
    min_breadth DOUBLE,
    min_ber DOUBLE,
    min_sylph_abundance DOUBLE,
    created_at DOUBLE NOT NULL,
    updated_at DOUBLE NOT NULL
);

CREATE TABLE IF NOT EXISTS matrix_store_samples (
    matrix_id VARCHAR NOT NULL,
    sample_id VARCHAR NOT NULL,
    added_at DOUBLE NOT NULL,
    PRIMARY KEY (matrix_id, sample_id)
);

CREATE TABLE IF NOT EXISTS matrix_compares (
    compare_id VARCHAR PRIMARY KEY,
    matrix_id VARCHAR NOT NULL,
    compare_db_file VARCHAR NOT NULL,
    calculate VARCHAR NOT NULL,
    created_at DOUBLE NOT NULL,
    updated_at DOUBLE NOT NULL
);

CREATE TABLE IF NOT EXISTS metatrawl_migrations (
    name VARCHAR PRIMARY KEY,
    applied_at DOUBLE NOT NULL
);
"""

DEFAULT_MEMORY_LIMIT_GB = 20.0
MEMORY_LIMIT_ENV_VAR = "METATRAWL_MEMORY_LIMIT_GB"
THREADS_ENV_VAR = "METATRAWL_THREADS"


@dataclass(frozen=True)
class ProfileBundle:
    """Completed output bundle for one SRA run."""

    run_id: str
    profile_file: Path
    genome_stats_file: Path
    sylph_abundance_file: Path
    gene_stats_file: Path | None = None


@dataclass(frozen=True)
class MatrixStore:
    """Registered matrix store metadata."""

    matrix_id: str
    genome: str
    matrix_file: Path
    profile_count: int
    storage_layout: str = "dense"


@dataclass(frozen=True)
class MatrixFilters:
    """Sample filters used when building a matrix store."""

    min_coverage: float | None = None
    min_breadth: float | None = None
    min_ber: float | None = None
    min_sylph_abundance: float | None = None


@dataclass(frozen=True)
class ProfileStorageConfig:
    """Immutable profile representation selected for one project database."""

    mode: str = allele_mask.PROFILE_STORAGE_FULL
    format_version: int = allele_mask.ALLELE_MASK_FORMAT_VERSION
    min_cov: int | None = None
    codec: str | None = None
    allele_bit_order: str | None = None


def connect(
    db_path: str | Path,
    *,
    memory_limit_gb: float | None = None,
    threads: int | None = None,
) -> duckdb.DuckDBPyConnection:
    """Open a DuckDB connection and ensure the MetaTrawl schema exists."""
    path = Path(db_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    conn = duckdb.connect(str(path))
    _apply_connection_settings(conn, memory_limit_gb=memory_limit_gb, threads=threads)
    init_schema(conn)
    return conn


def connect_read_only(
    db_path: str | Path,
    *,
    memory_limit_gb: float | None = None,
    threads: int | None = None,
) -> duckdb.DuckDBPyConnection:
    """Open a read-only DuckDB connection for concurrent matrix jobs."""
    conn = duckdb.connect(str(Path(db_path)), read_only=True)
    _apply_connection_settings(conn, memory_limit_gb=memory_limit_gb, threads=threads)
    return conn


def resolve_memory_limit_gb(memory_limit_gb: float | None = None) -> float | None:
    """Resolve the connection memory limit from the argument, env var, or default.

    A non-positive value disables the limit and restores DuckDB's own default.
    """
    if memory_limit_gb is None:
        raw = os.environ.get(MEMORY_LIMIT_ENV_VAR)
        if raw is not None and raw.strip():
            try:
                memory_limit_gb = float(raw)
            except ValueError as exc:
                raise ValueError(f"{MEMORY_LIMIT_ENV_VAR} must be a number, got: {raw!r}") from exc
        else:
            memory_limit_gb = DEFAULT_MEMORY_LIMIT_GB
    return memory_limit_gb if memory_limit_gb > 0 else None


def _resolve_threads(threads: int | None) -> int | None:
    if threads is not None:
        return max(1, int(threads))
    raw = os.environ.get(THREADS_ENV_VAR)
    if raw is None or not raw.strip():
        return None
    try:
        return max(1, int(raw))
    except ValueError as exc:
        raise ValueError(f"{THREADS_ENV_VAR} must be an integer, got: {raw!r}") from exc


def _apply_connection_settings(
    conn: duckdb.DuckDBPyConnection,
    *,
    memory_limit_gb: float | None,
    threads: int | None,
) -> None:
    """Bound DuckDB memory so long imports spill instead of exhausting the host."""
    resolved_memory = resolve_memory_limit_gb(memory_limit_gb)
    if resolved_memory is not None:
        conn.execute(f"SET memory_limit = {_sql_literal(f'{resolved_memory}GB')}")
    resolved_threads = _resolve_threads(threads)
    if resolved_threads is not None:
        conn.execute(f"SET threads = {resolved_threads}")


def init_schema(conn: duckdb.DuckDBPyConnection) -> None:
    """Create and lightly migrate registry tables."""
    had_legacy_matrix_profiles = _table_exists(conn, "matrix_store_profiles")
    conn.execute(SCHEMA_SQL)
    conn.execute("ALTER TABLE sra_runs ADD COLUMN IF NOT EXISTS status VARCHAR DEFAULT 'added'")
    conn.execute("ALTER TABLE profiles ADD COLUMN IF NOT EXISTS gene_stats_file VARCHAR")
    conn.execute(
        "ALTER TABLE profiles ADD COLUMN IF NOT EXISTS "
        "profile_storage_mode VARCHAR DEFAULT 'full'"
    )
    conn.execute("ALTER TABLE profiles ADD COLUMN IF NOT EXISTS profile_min_cov INTEGER")
    conn.execute("ALTER TABLE profile_positions ADD COLUMN IF NOT EXISTS ref_base_bitmask UTINYINT")
    conn.execute("ALTER TABLE genome_stats ADD COLUMN IF NOT EXISTS ref_ani DOUBLE")
    conn.execute("ALTER TABLE genome_stats ADD COLUMN IF NOT EXISTS coverage_median DOUBLE")
    conn.execute("ALTER TABLE genome_stats ADD COLUMN IF NOT EXISTS coverage_std DOUBLE")
    conn.execute("ALTER TABLE genome_stats ADD COLUMN IF NOT EXISTS genome_length BIGINT")
    conn.execute("ALTER TABLE genome_stats ADD COLUMN IF NOT EXISTS gap_mean DOUBLE")
    conn.execute("ALTER TABLE genome_stats ADD COLUMN IF NOT EXISTS gap_std DOUBLE")
    conn.execute('ALTER TABLE genome_stats ADD COLUMN IF NOT EXISTS "5x_cov_sites" BIGINT')
    conn.execute("ALTER TABLE genome_stats ADD COLUMN IF NOT EXISTS heterogeneity DOUBLE")
    conn.execute("ALTER TABLE genome_stats ADD COLUMN IF NOT EXISTS fug DOUBLE")
    conn.execute("ALTER TABLE genome_stats ADD COLUMN IF NOT EXISTS reads_mapped BIGINT")
    conn.execute("ALTER TABLE genome_stats ADD COLUMN IF NOT EXISTS conANI_reference DOUBLE")
    conn.execute("ALTER TABLE genome_stats ADD COLUMN IF NOT EXISTS SNS_count BIGINT")
    conn.execute("ALTER TABLE genome_stats ADD COLUMN IF NOT EXISTS SNV_count BIGINT")
    conn.execute("ALTER TABLE genome_stats ADD COLUMN IF NOT EXISTS presence VARCHAR")
    conn.execute("ALTER TABLE genome_stats ADD COLUMN IF NOT EXISTS genome_taxonomy VARCHAR")
    conn.execute("ALTER TABLE gene_stats ADD COLUMN IF NOT EXISTS ref_ani DOUBLE")
    conn.execute("ALTER TABLE gene_stats ADD COLUMN IF NOT EXISTS length BIGINT")
    conn.execute("ALTER TABLE matrix_stores ADD COLUMN IF NOT EXISTS storage_layout VARCHAR DEFAULT 'dense'")
    conn.execute("ALTER TABLE matrix_stores ADD COLUMN IF NOT EXISTS min_coverage DOUBLE")
    conn.execute("ALTER TABLE matrix_stores ADD COLUMN IF NOT EXISTS min_breadth DOUBLE")
    conn.execute("ALTER TABLE matrix_stores ADD COLUMN IF NOT EXISTS min_ber DOUBLE")
    conn.execute("ALTER TABLE matrix_stores ADD COLUMN IF NOT EXISTS min_sylph_abundance DOUBLE")
    conn.execute(
        """
        INSERT INTO profile_storage
          (id, mode, format_version, min_cov, codec, allele_bit_order)
        SELECT 1, 'full', ?, NULL, NULL, NULL
        WHERE NOT EXISTS (SELECT 1 FROM profile_storage WHERE id = 1)
        """,
        [allele_mask.ALLELE_MASK_FORMAT_VERSION],
    )
    # MetaTrawl never reads profile_positions.gene: the matrix export selects
    # chrom/pos/ACGT only. Dropping it is metadata-only and idempotent.
    conn.execute("ALTER TABLE profile_positions DROP COLUMN IF EXISTS gene")
    _migrate_profile_counts_to_uint16(conn)
    _run_once(conn, "normalize_sylph_genomes", _normalize_existing_sylph_genomes)
    if had_legacy_matrix_profiles:
        conn.execute(
            """
            INSERT OR IGNORE INTO matrix_store_samples
            SELECT matrix_id, run_id AS sample_id, added_at
            FROM matrix_store_profiles
            """
        )


def profile_storage_config(conn: duckdb.DuckDBPyConnection) -> ProfileStorageConfig:
    """Return the database profile-storage contract.

    Databases created by older MetaTrawl versions have no contract table and are
    therefore interpreted as full-profile databases.
    """
    if not _table_exists(conn, "profile_storage"):
        return ProfileStorageConfig()
    row = conn.execute(
        """
        SELECT mode, format_version, min_cov, codec, allele_bit_order
        FROM profile_storage
        WHERE id = 1
        """
    ).fetchone()
    if row is None:
        return ProfileStorageConfig()
    mode = str(row[0])
    if mode not in allele_mask.PROFILE_STORAGE_MODES:
        raise ValueError(f"Unsupported database profile storage mode: {mode}")
    config = ProfileStorageConfig(
        mode=mode,
        format_version=int(row[1]),
        min_cov=int(row[2]) if row[2] is not None else None,
        codec=str(row[3]) if row[3] is not None else None,
        allele_bit_order=str(row[4]) if row[4] is not None else None,
    )
    if mode == allele_mask.PROFILE_STORAGE_ALLELE_MASK:
        if config.format_version != allele_mask.ALLELE_MASK_FORMAT_VERSION:
            raise ValueError(
                "Unsupported allele-mask format version: "
                f"{config.format_version}; expected {allele_mask.ALLELE_MASK_FORMAT_VERSION}."
            )
        if config.min_cov is None or config.min_cov < 1:
            raise ValueError("Allele-mask database has no valid fixed min_cov.")
        if config.codec != allele_mask.ALLELE_MASK_CODEC:
            raise ValueError(f"Unsupported allele-mask codec: {config.codec}")
        if config.allele_bit_order != allele_mask.ALLELE_MASK_BIT_ORDER:
            raise ValueError(
                f"Unsupported allele-mask bit order: {config.allele_bit_order}"
            )
    return config


def configure_profile_storage(
    conn: duckdb.DuckDBPyConnection,
    *,
    mode: str,
    min_cov: int | None = None,
) -> ProfileStorageConfig:
    """Set profile storage before any samples have been imported."""
    if mode not in allele_mask.PROFILE_STORAGE_MODES:
        raise ValueError(
            "profile storage must be one of: "
            + ", ".join(allele_mask.PROFILE_STORAGE_MODES)
        )
    if mode == allele_mask.PROFILE_STORAGE_ALLELE_MASK:
        if min_cov is None or min_cov < 1:
            raise ValueError("Allele-mask profile storage requires --profile-min-cov >= 1.")
    elif min_cov is not None:
        raise ValueError("--profile-min-cov can only be used with allele-mask storage.")

    current = profile_storage_config(conn)
    requested = ProfileStorageConfig(
        mode=mode,
        format_version=allele_mask.ALLELE_MASK_FORMAT_VERSION,
        min_cov=min_cov,
        codec=(
            allele_mask.ALLELE_MASK_CODEC
            if mode == allele_mask.PROFILE_STORAGE_ALLELE_MASK
            else None
        ),
        allele_bit_order=(
            allele_mask.ALLELE_MASK_BIT_ORDER
            if mode == allele_mask.PROFILE_STORAGE_ALLELE_MASK
            else None
        ),
    )
    if current == requested:
        return current
    stored_profiles = int(
        conn.execute(
            """
            SELECT
              (SELECT count(*) FROM samples)
              + (SELECT count(*) FROM profile_positions)
              + (SELECT count(*) FROM allele_mask_profile_blocks)
            """
        ).fetchone()[0]
    )
    if stored_profiles:
        raise ValueError(
            f"Cannot change profile storage from {current.mode} to {mode} after "
            "samples have been imported. Create a new database instead."
        )
    conn.execute(
        """
        INSERT OR REPLACE INTO profile_storage
          (id, mode, format_version, min_cov, codec, allele_bit_order)
        VALUES (1, ?, ?, ?, ?, ?)
        """,
        [
            requested.mode,
            requested.format_version,
            requested.min_cov,
            requested.codec,
            requested.allele_bit_order,
        ],
    )
    return requested


def add_runs(conn: duckdb.DuckDBPyConnection, run_ids: list[str]) -> tuple[int, int]:
    """Add or reactivate SRA run IDs. Returns ``(added, reactivated)``."""
    now = time.time()
    added = 0
    reactivated = 0
    for run_id in _normalize_ids(run_ids):
        existing = conn.execute("SELECT deleted_at FROM sra_runs WHERE run_id = ?", [run_id]).fetchone()
        if existing is None:
            conn.execute(
                "INSERT INTO sra_runs (run_id, status, added_at, updated_at, deleted_at) VALUES (?, 'added', ?, ?, NULL)",
                [run_id, now, now],
            )
            added += 1
        elif existing[0] is not None:
            conn.execute(
                "UPDATE sra_runs SET deleted_at = NULL, status = 'added', updated_at = ? WHERE run_id = ?",
                [now, run_id],
            )
            reactivated += 1
    return added, reactivated


def delete_runs(conn: duckdb.DuckDBPyConnection, run_ids: list[str]) -> int:
    """Soft-delete active SRA runs and return the number changed."""
    now = time.time()
    changed = 0
    for run_id in _normalize_ids(run_ids):
        row = conn.execute(
            """
            UPDATE sra_runs
            SET deleted_at = ?, status = 'deleted', updated_at = ?
            WHERE run_id = ? AND deleted_at IS NULL
            RETURNING run_id
            """,
            [now, now, run_id],
        ).fetchone()
        if row is not None:
            changed += 1
    return changed


def list_runs(conn: duckdb.DuckDBPyConnection, include_deleted: bool = False) -> list[dict[str, object]]:
    """List registered SRA runs."""
    where = "" if include_deleted else "WHERE deleted_at IS NULL"
    return _rows_as_dicts(
        conn.execute(
            f"""
            SELECT run_id, status, added_at, updated_at, deleted_at
            FROM sra_runs
            {where}
            ORDER BY run_id
            """
        )
    )


def remaining_runs(conn: duckdb.DuckDBPyConnection) -> list[str]:
    """Return active SRA runs without a completed sample import."""
    rows = conn.execute(
        """
        SELECT r.run_id
        FROM sra_runs r
        LEFT JOIN samples s
          ON s.run_id = r.run_id AND s.status = 'complete'
        WHERE r.deleted_at IS NULL
          AND s.sample_id IS NULL
        ORDER BY r.run_id
        """
    ).fetchall()
    return [str(row[0]) for row in rows]


def mark_run_failed(conn: duckdb.DuckDBPyConnection, *, run_id: str, error: str) -> None:
    """Mark a run as failed while leaving it eligible for the next sync."""
    conn.execute(
        "UPDATE sra_runs SET status = 'failed', updated_at = ? WHERE run_id = ?",
        [time.time(), run_id],
    )


def import_profile_bundle(
    conn: duckdb.DuckDBPyConnection,
    bundle: ProfileBundle,
    *,
    add_run_if_missing: bool = False,
    cache_dir: Path | None = None,
    reference_cache: allele_mask.AlleleMaskWriteCache | None = None,
) -> None:
    """Import profile, stat, and abundance files into project tables."""
    import_profile_bundles(
        conn,
        [bundle],
        add_runs_if_missing=add_run_if_missing,
        cache_dir=cache_dir,
        reference_cache=reference_cache,
    )


def import_profile_bundles(
    conn: duckdb.DuckDBPyConnection,
    bundles: list[ProfileBundle],
    *,
    add_runs_if_missing: bool = False,
    cache_dir: Path | None = None,
    reference_cache: allele_mask.AlleleMaskWriteCache | None = None,
) -> None:
    """Import multiple bundles in one transaction using one DuckDB writer."""
    if not bundles:
        return
    for bundle in bundles:
        _validate_profile_bundle(
            conn,
            bundle,
            add_run_if_missing=add_runs_if_missing,
        )
    storage = profile_storage_config(conn)
    conn.execute("BEGIN TRANSACTION")
    try:
        for bundle in bundles:
            _import_profile_bundle_rows(
                conn,
                bundle,
                storage=storage,
                cache_dir=cache_dir,
                reference_cache=reference_cache,
            )
        conn.execute("COMMIT")
    except Exception:
        _rollback_quietly(conn)
        raise


def _validate_profile_bundle(
    conn: duckdb.DuckDBPyConnection,
    bundle: ProfileBundle,
    *,
    add_run_if_missing: bool,
) -> None:
    if conn.execute("SELECT 1 FROM sra_runs WHERE run_id = ?", [bundle.run_id]).fetchone() is None:
        if not add_run_if_missing:
            raise ValueError(f"Cannot import profile for unknown run_id: {bundle.run_id}")
        add_runs(conn, [bundle.run_id])
    _require_existing_file(bundle.profile_file, "profile_file")
    _require_existing_file(bundle.genome_stats_file, "genome_stats_file")
    _require_existing_file(bundle.sylph_abundance_file, "sylph_abundance_file")
    if bundle.gene_stats_file is not None:
        _require_existing_file(bundle.gene_stats_file, "gene_stats_file")


def _import_profile_bundle_rows(
    conn: duckdb.DuckDBPyConnection,
    bundle: ProfileBundle,
    *,
    storage: ProfileStorageConfig,
    cache_dir: Path | None,
    reference_cache: allele_mask.AlleleMaskWriteCache | None = None,
) -> None:
    sample_id = bundle.run_id
    now = time.time()
    existing = conn.execute(
        "SELECT created_at FROM samples WHERE sample_id = ?",
        [sample_id],
    ).fetchone()
    if existing is not None:
        _delete_existing_profile_rows(conn, sample_id=sample_id)

    if storage.mode == allele_mask.PROFILE_STORAGE_ALLELE_MASK:
        allele_mask.store_profile_parquet(
            conn,
            sample_id=sample_id,
            profile_file=bundle.profile_file,
            min_cov=int(storage.min_cov),
            cache_dir=cache_dir,
            reference_cache=reference_cache,
        )
    else:
        _insert_profile_positions(
            conn,
            sample_id=sample_id,
            profile_file=bundle.profile_file,
        )
    _insert_genome_stats(conn, sample_id=sample_id, stats_file=bundle.genome_stats_file)
    if bundle.gene_stats_file is not None:
        _insert_gene_stats(conn, sample_id=sample_id, stats_file=bundle.gene_stats_file)
    _insert_sylph_abundance(conn, sample_id=sample_id, abundance_file=bundle.sylph_abundance_file)

    created_at = float(existing[0]) if existing is not None else now
    conn.execute(
        "INSERT OR REPLACE INTO samples VALUES (?, ?, 'complete', ?, ?)",
        [sample_id, bundle.run_id, created_at, now],
    )
    profile_existing = conn.execute("SELECT created_at FROM profiles WHERE run_id = ?", [bundle.run_id]).fetchone()
    profile_created_at = float(profile_existing[0]) if profile_existing is not None else now
    conn.execute(
        """
        INSERT OR REPLACE INTO profiles
          (run_id, profile_file, genome_stats_file, gene_stats_file,
           sylph_abundance_file, profile_storage_mode, profile_min_cov,
           created_at, updated_at)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        [
            bundle.run_id,
            str(bundle.profile_file),
            str(bundle.genome_stats_file),
            str(bundle.gene_stats_file) if bundle.gene_stats_file is not None else None,
            str(bundle.sylph_abundance_file),
            storage.mode,
            storage.min_cov,
            profile_created_at,
            now,
        ],
    )
    conn.execute(
        "UPDATE sra_runs SET status = 'complete', updated_at = ? WHERE run_id = ?",
        [now, bundle.run_id],
    )


def _delete_existing_profile_rows(
    conn: duckdb.DuckDBPyConnection,
    *,
    sample_id: str,
) -> None:
    """Remove old payload rows only when explicitly replacing a sample."""
    conn.execute("DELETE FROM profile_positions WHERE sample_id = ?", [sample_id])
    conn.execute(
        "DELETE FROM allele_mask_profile_blocks WHERE sample_id = ?",
        [sample_id],
    )
    conn.execute("DELETE FROM genome_stats WHERE sample_id = ?", [sample_id])
    conn.execute("DELETE FROM gene_stats WHERE sample_id = ?", [sample_id])
    conn.execute("DELETE FROM sylph_abundance WHERE sample_id = ?", [sample_id])


def _rollback_quietly(conn: duckdb.DuckDBPyConnection) -> None:
    try:
        conn.execute("ROLLBACK")
    except Exception:
        pass


def add_profiles(
    conn: duckdb.DuckDBPyConnection,
    bundles: list[ProfileBundle],
    *,
    add_runs_if_missing: bool = False,
    cache_dir: Path | None = None,
) -> int:
    """Import completed profile bundles by SRA run ID."""
    for bundle in bundles:
        import_profile_bundle(
            conn,
            bundle,
            add_run_if_missing=add_runs_if_missing,
            cache_dir=cache_dir,
        )
    return len(bundles)


def list_profiles(conn: duckdb.DuckDBPyConnection) -> list[dict[str, object]]:
    """List completed samples and their source files."""
    return _rows_as_dicts(
        conn.execute(
            """
            SELECT run_id, profile_file, genome_stats_file, gene_stats_file,
                   sylph_abundance_file, profile_storage_mode, profile_min_cov,
                   created_at, updated_at
            FROM profiles
            ORDER BY run_id
            """
        )
    )


def eligible_sample_ids(conn: duckdb.DuckDBPyConnection, *, genome: str, filters: MatrixFilters) -> list[str]:
    """Return samples eligible for a matrix build under genome/stat thresholds."""
    conditions = ["s.status = 'complete'", "gs.genome = ?"]
    params: list[object] = [genome]
    if filters.min_coverage is not None:
        conditions.append("COALESCE(gs.coverage, 0) >= ?")
        params.append(filters.min_coverage)
    if filters.min_breadth is not None:
        conditions.append("COALESCE(gs.breadth, 0) >= ?")
        params.append(filters.min_breadth)
    if filters.min_ber is not None:
        conditions.append("COALESCE(gs.ber, 0) >= ?")
        params.append(filters.min_ber)
    if filters.min_sylph_abundance is not None:
        conditions.append(
            """
            EXISTS (
              SELECT 1
              FROM sylph_abundance sa
              WHERE sa.sample_id = s.sample_id
                AND (sa.genome = ? OR sa.accession = ?)
                AND COALESCE(sa.abundance, 0) >= ?
            )
            """
        )
        params.extend([genome, genome, filters.min_sylph_abundance])
    rows = conn.execute(
        f"""
        SELECT DISTINCT s.sample_id
        FROM samples s
        JOIN genome_stats gs USING (sample_id)
        WHERE {' AND '.join(conditions)}
        ORDER BY s.sample_id
        """,
        params,
    ).fetchall()
    return [str(row[0]) for row in rows]


def eligible_unmaterialized_sample_ids(
    conn: duckdb.DuckDBPyConnection,
    *,
    matrix_id: str | None = None,
    existing_sample_ids: list[str] | None = None,
    genome: str,
    filters: MatrixFilters,
) -> list[str]:
    """Return matrix-eligible samples not yet materialized into a matrix store."""
    eligible = eligible_sample_ids(conn, genome=genome, filters=filters)
    if not eligible:
        return []
    existing_sample_ids = existing_sample_ids or []
    if matrix_id is None:
        return [sample_id for sample_id in eligible if sample_id not in set(existing_sample_ids)]
    rows = conn.execute(
        """
        SELECT sample_id
        FROM unnest(?) AS candidates(sample_id)
        WHERE sample_id NOT IN (
            SELECT sample_id
            FROM matrix_store_samples
            WHERE matrix_id = ?
        )
        ORDER BY sample_id
        """,
        [eligible, matrix_id],
    ).fetchall()
    return [str(row[0]) for row in rows]


def genomes_with_unmaterialized_samples(
    conn: duckdb.DuckDBPyConnection,
    *,
    targets: list[tuple[str, str]],
    filters: MatrixFilters,
) -> list[str]:
    """Return target genomes holding eligible samples their matrix does not have.

    ``targets`` pairs each genome with the matrix ID that stores it. One grouped
    query replaces a per-genome eligibility scan, and lets the caller skip
    launching a job for a genome that has nothing new to materialize.
    """
    if not targets:
        return []
    conditions = ["s.status = 'complete'"]
    params: list[object] = []
    if filters.min_coverage is not None:
        conditions.append("COALESCE(gs.coverage, 0) >= ?")
        params.append(filters.min_coverage)
    if filters.min_breadth is not None:
        conditions.append("COALESCE(gs.breadth, 0) >= ?")
        params.append(filters.min_breadth)
    if filters.min_ber is not None:
        conditions.append("COALESCE(gs.ber, 0) >= ?")
        params.append(filters.min_ber)
    if filters.min_sylph_abundance is not None:
        conditions.append(
            """
            EXISTS (
              SELECT 1
              FROM sylph_abundance sa
              WHERE sa.sample_id = s.sample_id
                AND (sa.genome = t.genome OR sa.accession = t.genome)
                AND COALESCE(sa.abundance, 0) >= ?
            )
            """
        )
        params.append(filters.min_sylph_abundance)
    frame = pl.DataFrame(
        {
            "genome": [str(genome) for genome, _matrix_id in targets],
            "matrix_id": [str(matrix_id) for _genome, matrix_id in targets],
        }
    )
    conn.register("_metatrawl_matrix_targets", frame)
    try:
        rows = conn.execute(
            f"""
            SELECT DISTINCT t.genome
            FROM _metatrawl_matrix_targets t
            JOIN genome_stats gs ON gs.genome = t.genome
            JOIN samples s ON s.sample_id = gs.sample_id
            WHERE {' AND '.join(conditions)}
              AND NOT EXISTS (
                SELECT 1
                FROM matrix_store_samples m
                WHERE m.matrix_id = t.matrix_id AND m.sample_id = s.sample_id
              )
            ORDER BY t.genome
            """,
            params,
        ).fetchall()
    finally:
        conn.unregister("_metatrawl_matrix_targets")
    return [str(row[0]) for row in rows]


def completed_sample_ids(conn: duckdb.DuckDBPyConnection) -> list[str]:
    """Return complete sample IDs."""
    rows = conn.execute(
        """
        SELECT sample_id
        FROM samples
        WHERE status = 'complete'
        ORDER BY sample_id
        """
    ).fetchall()
    return [str(row[0]) for row in rows]


def genomes_with_complete_samples(conn: duckdb.DuckDBPyConnection) -> list[str]:
    """Return genomes represented by at least one complete sample."""
    rows = conn.execute(
        """
        SELECT DISTINCT gs.genome
        FROM genome_stats gs
        JOIN samples s USING (sample_id)
        WHERE s.status = 'complete'
        ORDER BY gs.genome
        """
    ).fetchall()
    return [str(row[0]) for row in rows]


def matrix_store_filters(conn: duckdb.DuckDBPyConnection, matrix_id: str) -> MatrixFilters:
    """Return the filters originally used to build a matrix store."""
    row = conn.execute(
        """
        SELECT min_coverage, min_breadth, min_ber, min_sylph_abundance
        FROM matrix_stores
        WHERE matrix_id = ?
        """,
        [matrix_id],
    ).fetchone()
    if row is None:
        raise ValueError(f"Unknown matrix ID: {matrix_id}")
    return MatrixFilters(
        min_coverage=row[0],
        min_breadth=row[1],
        min_ber=row[2],
        min_sylph_abundance=row[3],
    )



def export_profile_parquets(
    conn: duckdb.DuckDBPyConnection,
    *,
    sample_ids: list[str],
    output_dir: Path,
    genome: str | None = None,
    memory_limit_gb: float | None = None,
    export_batch_mb: float = 128.0,
    duckdb_threads: int = 1,
    progress_callback: ExportProgressCallback | None = None,
) -> list[Path]:
    """Export selected samples from DuckDB into temporary ZipStrain profile parquets."""
    storage = profile_storage_config(conn)
    if storage.mode == allele_mask.PROFILE_STORAGE_ALLELE_MASK:
        raise ValueError(
            "Allele-mask databases do not contain reconstructable A/C/G/T counts. "
            "Build a bitmask matrix directly instead."
        )
    output_dir.mkdir(parents=True, exist_ok=True)
    _configure_duckdb_export(conn, output_dir=output_dir, memory_limit_gb=memory_limit_gb, threads=duckdb_threads)
    row_group_size = _export_row_group_size(export_batch_mb)
    paths: list[Path] = []
    total = len(sample_ids)
    if progress_callback is not None:
        progress_callback({"phase": "start", "completed": 0, "total": total, "genome": genome or "all"})
    for index, sample_id in enumerate(sample_ids, start=1):
        output_file = output_dir / f"{sample_id}.parquet"
        conditions = [f"sample_id = {_sql_literal(sample_id)}"]
        count_params: list[object] = [sample_id]
        count_conditions = ["sample_id = ?"]
        if genome is not None and genome != "all":
            conditions.append(f"genome = {_sql_literal(genome)}")
            count_conditions.append("genome = ?")
            count_params.append(genome)
        count_where_sql = " AND ".join(count_conditions)
        row_count = conn.execute(
            f"SELECT count(*) FROM profile_positions WHERE {count_where_sql}",
            count_params,
        ).fetchone()[0]
        if row_count:
            where_sql = " AND ".join(conditions)
            conn.execute(
                f"""
                COPY (
                    SELECT chrom, genome, pos, COALESCE(gene, 'NA') AS gene, A, T, C, G
                    FROM profile_positions
                    WHERE {where_sql}
                ) TO {_sql_literal(output_file)} (FORMAT PARQUET, ROW_GROUP_SIZE {row_group_size})
                """
            )
            paths.append(output_file)
        if progress_callback is not None:
            progress_callback(
                {
                    "phase": "advance",
                    "completed": index,
                    "total": total,
                    "sample_name": sample_id,
                    "genome": genome or "all",
                }
            )
    if progress_callback is not None:
        progress_callback({"phase": "done", "completed": total, "total": total, "genome": genome or "all"})
    return paths


def _configure_duckdb_export(
    conn: duckdb.DuckDBPyConnection,
    *,
    output_dir: Path,
    memory_limit_gb: float | None,
    threads: int,
) -> None:
    # Matrix export can otherwise use large parallel parquet buffers per sample.
    conn.execute(f"SET threads = {max(1, int(threads))}")
    temp_dir = output_dir / ".duckdb_tmp"
    temp_dir.mkdir(parents=True, exist_ok=True)
    conn.execute(f"SET temp_directory = {_sql_literal(temp_dir)}")
    if memory_limit_gb is not None:
        conn.execute(f"SET memory_limit = {_sql_literal(f'{memory_limit_gb}GB')}")


def _export_row_group_size(export_batch_mb: float) -> int:
    estimated_row_width = 128
    target_rows = int(max(1.0, export_batch_mb) * 1024 * 1024 / estimated_row_width)
    return max(10_000, min(250_000, target_rows))


def _sql_literal(value: object) -> str:
    return "'" + str(value).replace("'", "''") + "'"


def register_cache_genome(conn: duckdb.DuckDBPyConnection, *, accession: str, genome_fasta: Path, gene_fasta: Path) -> None:
    """Record a prepared cache genome."""
    conn.execute(
        "INSERT OR REPLACE INTO cache_genomes VALUES (?, ?, ?, 'ready', ?)",
        [accession, str(genome_fasta), str(gene_fasta), time.time()],
    )


def register_matrix_store(
    conn: duckdb.DuckDBPyConnection,
    *,
    matrix_id: str,
    genome: str,
    matrix_file: Path,
    profile_count: int,
    storage_layout: str = "dense",
    sample_ids: list[str] | None = None,
    filters: MatrixFilters | None = None,
    overwrite: bool = False,
) -> MatrixStore:
    """Record matrix metadata for legacy listing; the HDF5 file is authoritative."""
    existing = get_matrix_store(conn, matrix_id)
    now = time.time()
    created_at = now
    if existing is not None:
        created_at = float(conn.execute("SELECT created_at FROM matrix_stores WHERE matrix_id = ?", [matrix_id]).fetchone()[0])
    filters = filters or MatrixFilters()
    conn.execute(
        """
        INSERT OR REPLACE INTO matrix_stores
          (matrix_id, genome, matrix_file, profile_count, storage_layout,
           min_coverage, min_breadth, min_ber, min_sylph_abundance, created_at, updated_at)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        [
            matrix_id,
            genome,
            str(matrix_file),
            profile_count,
            storage_layout,
            filters.min_coverage,
            filters.min_breadth,
            filters.min_ber,
            filters.min_sylph_abundance,
            created_at,
            now,
        ],
    )
    conn.execute("DELETE FROM matrix_store_samples WHERE matrix_id = ?", [matrix_id])
    if sample_ids:
        add_matrix_store_samples(conn, matrix_id=matrix_id, sample_ids=sample_ids)
    return MatrixStore(matrix_id=matrix_id, genome=genome, matrix_file=matrix_file, profile_count=profile_count, storage_layout=storage_layout)


def get_matrix_store(conn: duckdb.DuckDBPyConnection, matrix_id: str) -> MatrixStore | None:
    """Return a registered matrix store by ID, if present."""
    row = conn.execute(
        "SELECT matrix_id, genome, matrix_file, profile_count, storage_layout FROM matrix_stores WHERE matrix_id = ?",
        [matrix_id],
    ).fetchone()
    if row is None:
        return None
    return MatrixStore(row[0], row[1], Path(row[2]), int(row[3]), row[4])


def get_matrix_store_by_file(conn: duckdb.DuckDBPyConnection, matrix_file: Path) -> MatrixStore | None:
    """Return a registered matrix store by matrix file path, if present."""
    requested = str(Path(matrix_file))
    requested_resolved = str(Path(matrix_file).expanduser().resolve())
    rows = conn.execute(
        """
        SELECT matrix_id, genome, matrix_file, profile_count, storage_layout
        FROM matrix_stores
        WHERE matrix_file = ?
           OR matrix_file = ?
        """,
        [requested, requested_resolved],
    ).fetchall()
    if not rows:
        all_rows = conn.execute(
            "SELECT matrix_id, genome, matrix_file, profile_count, storage_layout FROM matrix_stores"
        ).fetchall()
        matches = [
            row
            for row in all_rows
            if Path(row[2]).expanduser().resolve() == Path(matrix_file).expanduser().resolve()
        ]
        rows = matches
    if not rows:
        return None
    if len(rows) > 1:
        raise ValueError(f"Multiple matrix registry rows point to matrix file: {matrix_file}")
    row = rows[0]
    return MatrixStore(row[0], row[1], Path(row[2]), int(row[3]), row[4])


def add_matrix_store_samples(conn: duckdb.DuckDBPyConnection, *, matrix_id: str, sample_ids: list[str]) -> int:
    """Record samples materialized into a matrix store.

    Registering a matrix with many samples is one bulk insert: row-at-a-time
    inserts cost hundreds of microseconds each, which a 100k-sample matrix turns
    into minutes.
    """
    normalized = _normalize_ids(sample_ids)
    if not normalized:
        return 0
    now = time.time()
    before = int(
        conn.execute(
            "SELECT count(*) FROM matrix_store_samples WHERE matrix_id = ?",
            [matrix_id],
        ).fetchone()[0]
    )
    frame = pl.DataFrame(
        {
            "matrix_id": [matrix_id] * len(normalized),
            "sample_id": normalized,
            "added_at": [now] * len(normalized),
        }
    )
    conn.register("_metatrawl_matrix_store_samples", frame)
    try:
        conn.execute(
            """
            INSERT OR IGNORE INTO matrix_store_samples (matrix_id, sample_id, added_at)
            SELECT matrix_id, sample_id, added_at FROM _metatrawl_matrix_store_samples
            """
        )
    finally:
        conn.unregister("_metatrawl_matrix_store_samples")
    after = int(
        conn.execute(
            "SELECT count(*) FROM matrix_store_samples WHERE matrix_id = ?",
            [matrix_id],
        ).fetchone()[0]
    )
    return after - before


def unmaterialized_sample_ids(conn: duckdb.DuckDBPyConnection, matrix_id: str) -> list[str]:
    """Return completed samples not yet materialized into a matrix."""
    rows = conn.execute(
        """
        SELECT s.sample_id
        FROM samples s
        LEFT JOIN matrix_store_samples mss
          ON mss.matrix_id = ? AND mss.sample_id = s.sample_id
        WHERE s.status = 'complete'
          AND mss.sample_id IS NULL
        ORDER BY s.sample_id
        """,
        [matrix_id],
    ).fetchall()
    return [str(row[0]) for row in rows]


def update_matrix_profile_count(conn: duckdb.DuckDBPyConnection, *, matrix_id: str) -> int:
    """Sync a matrix store's profile count from its materialized sample mapping."""
    count = int(conn.execute("SELECT count(*) FROM matrix_store_samples WHERE matrix_id = ?", [matrix_id]).fetchone()[0])
    conn.execute("UPDATE matrix_stores SET profile_count = ?, updated_at = ? WHERE matrix_id = ?", [count, time.time(), matrix_id])
    return count


def register_matrix_compare(
    conn: duckdb.DuckDBPyConnection,
    *,
    compare_id: str,
    matrix_id: str,
    compare_db_file: Path,
    calculate: str,
) -> None:
    """Record a matrix compare output in the registry."""
    now = time.time()
    existing = conn.execute("SELECT created_at FROM matrix_compares WHERE compare_id = ?", [compare_id]).fetchone()
    created_at = float(existing[0]) if existing is not None else now
    conn.execute(
        "INSERT OR REPLACE INTO matrix_compares VALUES (?, ?, ?, ?, ?, ?)",
        [compare_id, matrix_id, str(compare_db_file), calculate, created_at, now],
    )


def registry_status(conn: duckdb.DuckDBPyConnection) -> dict[str, object]:
    """Return high-level registry counts."""
    active_runs = conn.execute("SELECT count(*) FROM sra_runs WHERE deleted_at IS NULL").fetchone()[0]
    deleted_runs = conn.execute("SELECT count(*) FROM sra_runs WHERE deleted_at IS NOT NULL").fetchone()[0]
    samples = conn.execute("SELECT count(*) FROM samples WHERE status = 'complete'").fetchone()[0]
    profile_rows = conn.execute("SELECT count(*) FROM profile_positions").fetchone()[0]
    allele_mask_blocks = conn.execute(
        "SELECT count(*) FROM allele_mask_profile_blocks"
    ).fetchone()[0]
    allele_mask_positions = conn.execute(
        "SELECT COALESCE(sum(covered_positions), 0) FROM allele_mask_profile_blocks"
    ).fetchone()[0]
    matrices = conn.execute("SELECT count(*) FROM matrix_stores").fetchone()[0]
    compares = conn.execute("SELECT count(*) FROM matrix_compares").fetchone()[0]
    remaining = len(remaining_runs(conn))
    storage = profile_storage_config(conn)
    return {
        "active_runs": int(active_runs),
        "deleted_runs": int(deleted_runs),
        "complete_samples": int(samples),
        "profile_rows": int(profile_rows),
        "profile_storage": storage.mode,
        "profile_min_cov": storage.min_cov if storage.min_cov is not None else "NA",
        "allele_mask_blocks": int(allele_mask_blocks),
        "allele_mask_positions": int(allele_mask_positions),
        "remaining_profiles": remaining,
        "matrix_stores": int(matrices),
        "matrix_compares": int(compares),
    }


def _insert_profile_positions(conn: duckdb.DuckDBPyConnection, *, sample_id: str, profile_file: Path) -> None:
    if profile_file.suffix.lower() != ".parquet":
        _insert_profile_positions_via_polars(conn, sample_id=sample_id, profile_file=profile_file)
        return

    columns = set(_duckdb_parquet_schema(conn, profile_file))
    required = {"chrom", "pos", "genome", "A", "C", "G", "T"}
    missing = required - columns
    if missing:
        raise ValueError(f"profile_file missing required columns: {', '.join(sorted(missing))}")
    ref_expr = "CAST(ref_base_bitmask AS UTINYINT)" if "ref_base_bitmask" in columns else "CAST(NULL AS UTINYINT)"
    conn.execute(
        f"""
        INSERT INTO profile_positions
          (sample_id, chrom, pos, genome, A, C, G, T, ref_base_bitmask)
        SELECT
          ? AS sample_id,
          CAST(chrom AS VARCHAR) AS chrom,
          CAST(pos AS BIGINT) AS pos,
          CAST(genome AS VARCHAR) AS genome,
          CAST(A AS USMALLINT) AS A,
          CAST(C AS USMALLINT) AS C,
          CAST(G AS USMALLINT) AS G,
          CAST(T AS USMALLINT) AS T,
          {ref_expr} AS ref_base_bitmask
        FROM read_parquet(?)
        """,
        [sample_id, str(profile_file)],
    )


def _insert_profile_positions_via_polars(conn: duckdb.DuckDBPyConnection, *, sample_id: str, profile_file: Path) -> None:
    df = pl.read_parquet(profile_file)
    required = {"chrom", "pos", "genome", "A", "C", "G", "T"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"profile_file missing required columns: {', '.join(sorted(missing))}")
    df = df.select(
        pl.lit(sample_id).alias("sample_id"),
        pl.col("chrom").cast(pl.Utf8),
        pl.col("pos").cast(pl.Int64),
        pl.col("genome").cast(pl.Utf8),
        pl.col("A").cast(pl.UInt16),
        pl.col("C").cast(pl.UInt16),
        pl.col("G").cast(pl.UInt16),
        pl.col("T").cast(pl.UInt16),
        (
            pl.col("ref_base_bitmask").cast(pl.UInt8)
            if "ref_base_bitmask" in df.columns
            else pl.lit(None, dtype=pl.UInt8)
        ).alias("ref_base_bitmask"),
    )
    conn.register("_metatrawl_profile_positions", df)
    try:
        conn.execute(
            """INSERT INTO profile_positions
               (sample_id, chrom, pos, genome, A, C, G, T, ref_base_bitmask)
               SELECT sample_id, chrom, pos, genome, A, C, G, T, ref_base_bitmask
               FROM _metatrawl_profile_positions"""
        )
    finally:
        conn.unregister("_metatrawl_profile_positions")


def _duckdb_parquet_schema(conn: duckdb.DuckDBPyConnection, path: Path) -> dict[str, str]:
    """Return the parquet column names mapped to their DuckDB type names."""
    rows = conn.execute("DESCRIBE SELECT * FROM read_parquet(?)", [str(path)]).fetchall()
    return {str(row[0]): str(row[1]).upper() for row in rows}


def _sql_identifier(name: str) -> str:
    return '"' + str(name).replace('"', '""') + '"'


_FLOAT_PARQUET_TYPES = {"FLOAT", "DOUBLE", "REAL", "DECIMAL"}


def _is_float_type(type_name: str) -> bool:
    return any(type_name.startswith(candidate) for candidate in _FLOAT_PARQUET_TYPES)


def _optional_sql_expr(schema: dict[str, str], names: list[str], sql_type: str) -> str:
    """Mirror `_optional_*_expr`, but as SQL over `read_parquet`.

    Polars casts float -> integer by truncating toward zero while DuckDB rounds,
    so float sources are truncated explicitly to keep imported values identical.
    """
    for name in names:
        if name in schema:
            column = _sql_identifier(name)
            if sql_type == "BIGINT" and _is_float_type(schema[name]):
                return f"CAST(TRUNC({column}) AS BIGINT)"
            return f"CAST({column} AS {sql_type})"
    return f"CAST(NULL AS {sql_type})"


GENOME_COLUMN_CANDIDATES = ["genome", "genome_name", "reference", "accession"]

# (target column, source column candidates, SQL type). Shared by the parquet and
# the CSV/TSV import paths so the two cannot drift apart.
_GENOME_STATS_FIELDS: tuple[tuple[str, list[str], str], ...] = (
    ("coverage", ["coverage", "cov", "mean_coverage"], "DOUBLE"),
    ("breadth", ["breadth", "breadth_coverage", "breadth_cov"], "DOUBLE"),
    ("ber", ["ber", "BER"], "DOUBLE"),
    ("ref_ani", ["ref_ani", "reference_ani"], "DOUBLE"),
    ("coverage_median", ["coverage_median"], "DOUBLE"),
    ("coverage_std", ["coverage_std"], "DOUBLE"),
    ("genome_length", ["genome_length", "length"], "BIGINT"),
    ("gap_mean", ["gap_mean"], "DOUBLE"),
    ("gap_std", ["gap_std"], "DOUBLE"),
    ("5x_cov_sites", ["5x_cov_sites"], "BIGINT"),
    ("heterogeneity", ["heterogeneity"], "DOUBLE"),
    ("fug", ["fug", "FUG"], "DOUBLE"),
    ("reads_mapped", ["reads_mapped"], "BIGINT"),
    ("conANI_reference", ["conANI_reference", "conani_reference"], "DOUBLE"),
    ("SNS_count", ["SNS_count", "sns_count"], "BIGINT"),
    ("SNV_count", ["SNV_count", "snv_count"], "BIGINT"),
    ("presence", ["presence"], "VARCHAR"),
    ("genome_taxonomy", ["genome_taxonomy", "taxonomy"], "VARCHAR"),
)

_GENE_STATS_FIELDS: tuple[tuple[str, list[str], str], ...] = (
    ("coverage", ["coverage", "cov", "mean_coverage"], "DOUBLE"),
    ("breadth", ["breadth", "breadth_coverage", "breadth_cov"], "DOUBLE"),
    ("ber", ["ber", "BER"], "DOUBLE"),
    ("ref_ani", ["ref_ani", "reference_ani"], "DOUBLE"),
    ("length", ["length", "gene_length"], "BIGINT"),
)


def _polars_optional_expr(df: pl.DataFrame, names: list[str], sql_type: str) -> pl.Expr:
    if sql_type == "DOUBLE":
        return _optional_float_expr(df, names)
    if sql_type == "BIGINT":
        return _optional_int_expr(df, names)
    return _optional_string_expr(df, names)


def _insert_genome_stats(conn: duckdb.DuckDBPyConnection, *, sample_id: str, stats_file: Path) -> None:
    target_columns = ["sample_id", "genome"] + [field[0] for field in _GENOME_STATS_FIELDS]
    if stats_file.suffix.lower() == ".parquet":
        schema = _duckdb_parquet_schema(conn, stats_file)
        genome_col = _first_existing_name(schema, GENOME_COLUMN_CANDIDATES)
        if genome_col is None:
            raise ValueError("genome_stats_file missing a genome column")
        selects = [
            "? AS sample_id",
            f"CAST({_sql_identifier(genome_col)} AS VARCHAR) AS genome",
        ] + [
            f"{_optional_sql_expr(schema, names, sql_type)} AS {_sql_identifier(column)}"
            for column, names, sql_type in _GENOME_STATS_FIELDS
        ]
        conn.execute(
            f"""INSERT INTO genome_stats ({', '.join(_sql_identifier(c) for c in target_columns)})
                SELECT {', '.join(selects)}
                FROM read_parquet(?)""",
            [sample_id, str(stats_file)],
        )
        return

    df = _read_table(stats_file)
    genome_col = _first_existing(df, GENOME_COLUMN_CANDIDATES)
    if genome_col is None:
        raise ValueError("genome_stats_file missing a genome column")
    df = df.select(
        pl.lit(sample_id).alias("sample_id"),
        pl.col(genome_col).cast(pl.Utf8).alias("genome"),
        *[
            _polars_optional_expr(df, names, sql_type).alias(column)
            for column, names, sql_type in _GENOME_STATS_FIELDS
        ],
    )
    _insert_registered_frame(conn, "_metatrawl_genome_stats", df, table="genome_stats", columns=target_columns)


def _insert_gene_stats(conn: duckdb.DuckDBPyConnection, *, sample_id: str, stats_file: Path) -> None:
    target_columns = ["sample_id", "genome", "gene"] + [field[0] for field in _GENE_STATS_FIELDS]
    if stats_file.suffix.lower() == ".parquet":
        schema = _duckdb_parquet_schema(conn, stats_file)
        gene_col = _first_existing_name(schema, ["gene", "gene_id"])
        if gene_col is None:
            raise ValueError("gene_stats_file missing a gene column")
        genome_col = _first_existing_name(schema, GENOME_COLUMN_CANDIDATES)
        selects = [
            "? AS sample_id",
            (
                f"CAST({_sql_identifier(genome_col)} AS VARCHAR) AS genome"
                if genome_col is not None
                else "CAST(NULL AS VARCHAR) AS genome"
            ),
            f"CAST({_sql_identifier(gene_col)} AS VARCHAR) AS gene",
        ] + [
            f"{_optional_sql_expr(schema, names, sql_type)} AS {_sql_identifier(column)}"
            for column, names, sql_type in _GENE_STATS_FIELDS
        ]
        conn.execute(
            f"""INSERT INTO gene_stats ({', '.join(_sql_identifier(c) for c in target_columns)})
                SELECT {', '.join(selects)}
                FROM read_parquet(?)""",
            [sample_id, str(stats_file)],
        )
        return

    df = _read_table(stats_file)
    gene_col = _first_existing(df, ["gene", "gene_id"])
    if gene_col is None:
        raise ValueError("gene_stats_file missing a gene column")
    genome_col = _first_existing(df, GENOME_COLUMN_CANDIDATES)
    df = df.select(
        pl.lit(sample_id).alias("sample_id"),
        (pl.col(genome_col).cast(pl.Utf8) if genome_col else pl.lit(None, dtype=pl.Utf8)).alias("genome"),
        pl.col(gene_col).cast(pl.Utf8).alias("gene"),
        *[
            _polars_optional_expr(df, names, sql_type).alias(column)
            for column, names, sql_type in _GENE_STATS_FIELDS
        ],
    )
    _insert_registered_frame(conn, "_metatrawl_gene_stats", df, table="gene_stats", columns=target_columns)


def _insert_registered_frame(
    conn: duckdb.DuckDBPyConnection,
    view_name: str,
    df: pl.DataFrame,
    *,
    table: str,
    columns: list[str],
) -> None:
    column_sql = ", ".join(_sql_identifier(column) for column in columns)
    conn.register(view_name, df)
    try:
        conn.execute(f"INSERT INTO {table} ({column_sql}) SELECT {column_sql} FROM {view_name}")
    finally:
        conn.unregister(view_name)


def _insert_sylph_abundance(conn: duckdb.DuckDBPyConnection, *, sample_id: str, abundance_file: Path) -> None:
    df = _read_table(abundance_file)
    semantic_genome_col = _first_existing(df, ["genome", "genome_name"])
    genome_col = semantic_genome_col or _first_existing(df, ["reference", "Genome_file", "name"])
    accession_col = _first_existing(df, ["accession", "genome", "genome_name", "reference", "Genome_file", "name"])
    abundance_col = _first_existing(df, ["abundance", "relative_abundance", "Taxonomic_abundance", "ANI"])
    if genome_col is None and accession_col is None:
        raise ValueError("sylph_abundance_file missing a genome/accession column")
    source_col = accession_col or genome_col
    canonical_accession = (
        pl.col(source_col)
        .cast(pl.Utf8)
        .str.extract(ACCESSION_PATTERN.pattern, 1)
        .str.to_uppercase()
    )
    genome_expr = (
        pl.col(semantic_genome_col).cast(pl.Utf8)
        if semantic_genome_col is not None
        else canonical_accession
    )
    df = df.select(
        pl.lit(sample_id).alias("sample_id"),
        genome_expr.alias("genome"),
        canonical_accession.alias("accession"),
        (pl.col(abundance_col).cast(pl.Float64) if abundance_col else pl.lit(None, dtype=pl.Float64)).alias("abundance"),
    )
    if df["accession"].null_count():
        bad_values = (
            _read_table(abundance_file)
            .filter(canonical_accession.is_null())
            .get_column(source_col)
            .cast(pl.Utf8)
            .unique()
            .head(5)
            .to_list()
        )
        raise ValueError(
            "sylph_abundance_file contains genome values without a recognizable "
            f"GCF/GCA accession: {', '.join(bad_values)}"
        )
    conn.register("_metatrawl_sylph_abundance", df)
    try:
        conn.execute("INSERT INTO sylph_abundance SELECT * FROM _metatrawl_sylph_abundance")
    finally:
        conn.unregister("_metatrawl_sylph_abundance")


def _read_table(path: Path) -> pl.DataFrame:
    suffix = path.suffix.lower()
    if suffix == ".parquet":
        return pl.read_parquet(path)
    if suffix in {".tsv", ".tab"}:
        return pl.read_csv(path, separator="\t")
    return pl.read_csv(path)


def _optional_float_expr(df: pl.DataFrame, names: list[str]) -> pl.Expr:
    column = _first_existing(df, names)
    if column is None:
        return pl.lit(None, dtype=pl.Float64)
    return pl.col(column).cast(pl.Float64)


def _optional_int_expr(df: pl.DataFrame, names: list[str]) -> pl.Expr:
    column = _first_existing(df, names)
    if column is None:
        return pl.lit(None, dtype=pl.Int64)
    return pl.col(column).cast(pl.Int64)


def _optional_string_expr(df: pl.DataFrame, names: list[str]) -> pl.Expr:
    column = _first_existing(df, names)
    if column is None:
        return pl.lit(None, dtype=pl.Utf8)
    return pl.col(column).cast(pl.Utf8)


def _first_existing(df: pl.DataFrame, names: list[str]) -> str | None:
    return _first_existing_name(df.columns, names)


def _first_existing_name(available: Iterable[str], names: list[str]) -> str | None:
    columns = set(available)
    for name in names:
        if name in columns:
            return name
    return None


def _normalize_ids(values: list[str]) -> list[str]:
    seen: set[str] = set()
    normalized: list[str] = []
    for value in values:
        clean = str(value).strip()
        if clean and clean not in seen:
            normalized.append(clean)
            seen.add(clean)
    return normalized


def _require_existing_file(path: Path, column: str) -> None:
    if not path.exists():
        raise FileNotFoundError(f"{column} does not exist: {path}")


def _table_exists(conn: duckdb.DuckDBPyConnection, table_name: str) -> bool:
    row = conn.execute(
        """
        SELECT 1
        FROM information_schema.tables
        WHERE table_schema = current_schema()
          AND table_name = ?
        """,
        [table_name],
    ).fetchone()
    return row is not None


def _migrate_profile_counts_to_uint16(conn: duckdb.DuckDBPyConnection) -> None:
    """Convert legacy floating-point profile counts to unsigned 16-bit integers."""
    column_types = {
        str(name): str(data_type).upper()
        for name, data_type in conn.execute(
            """
            SELECT column_name, data_type
            FROM information_schema.columns
            WHERE table_schema = current_schema()
              AND table_name = 'profile_positions'
              AND column_name IN ('A', 'C', 'G', 'T')
            """
        ).fetchall()
    }
    for column in ("A", "C", "G", "T"):
        if column_types.get(column) != "USMALLINT":
            conn.execute(
                f'ALTER TABLE profile_positions ALTER COLUMN "{column}" '
                f'TYPE USMALLINT USING "{column}"::USMALLINT'
            )


def _run_once(
    conn: duckdb.DuckDBPyConnection,
    name: str,
    migration: Callable[[duckdb.DuckDBPyConnection], None],
) -> bool:
    """Run a one-time migration and record it so later connections skip the scan."""
    applied = conn.execute(
        "SELECT 1 FROM metatrawl_migrations WHERE name = ?",
        [name],
    ).fetchone()
    if applied is not None:
        return False
    migration(conn)
    conn.execute(
        "INSERT OR REPLACE INTO metatrawl_migrations (name, applied_at) VALUES (?, ?)",
        [name, time.time()],
    )
    return True


def _normalize_existing_sylph_genomes(conn: duckdb.DuckDBPyConnection) -> None:
    """Replace legacy Sylph database paths with canonical assembly accessions.

    Rows written by `_insert_sylph_abundance` are already canonical, so this only
    needs to run once per database.
    """
    pattern = ACCESSION_PATTERN.pattern.replace("'", "''")
    conn.execute(
        f"""
        UPDATE sylph_abundance
        SET genome = upper(regexp_extract(genome, '{pattern}', 1))
        WHERE regexp_extract(genome, '{pattern}', 1) <> ''
          AND genome <> upper(regexp_extract(genome, '{pattern}', 1))
        """
    )
    conn.execute(
        f"""
        UPDATE sylph_abundance
        SET accession = upper(regexp_extract(accession, '{pattern}', 1))
        WHERE regexp_extract(accession, '{pattern}', 1) <> ''
          AND accession <> upper(regexp_extract(accession, '{pattern}', 1))
        """
    )


def _rows_as_dicts(result) -> list[dict[str, object]]:
    columns = [desc[0] for desc in result.description]
    return [dict(zip(columns, row)) for row in result.fetchall()]

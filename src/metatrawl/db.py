"""DuckDB project database for MetaTrawl."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re
import time
from typing import Callable

import duckdb
import polars as pl


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
    created_at DOUBLE NOT NULL,
    updated_at DOUBLE NOT NULL
);

CREATE TABLE IF NOT EXISTS profile_positions (
    sample_id VARCHAR NOT NULL,
    chrom VARCHAR NOT NULL,
    pos BIGINT NOT NULL,
    genome VARCHAR NOT NULL,
    gene VARCHAR,
    A USMALLINT NOT NULL,
    C USMALLINT NOT NULL,
    G USMALLINT NOT NULL,
    T USMALLINT NOT NULL,
    ref_base_bitmask UTINYINT
);

CREATE TABLE IF NOT EXISTS genome_stats (
    sample_id VARCHAR NOT NULL,
    genome VARCHAR NOT NULL,
    coverage DOUBLE,
    breadth DOUBLE,
    ber DOUBLE,
    ref_ani DOUBLE
);

CREATE TABLE IF NOT EXISTS gene_stats (
    sample_id VARCHAR NOT NULL,
    genome VARCHAR,
    gene VARCHAR NOT NULL,
    coverage DOUBLE,
    breadth DOUBLE,
    ber DOUBLE,
    ref_ani DOUBLE
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
"""


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


def connect(db_path: str | Path) -> duckdb.DuckDBPyConnection:
    """Open a DuckDB connection and ensure the MetaTrawl schema exists."""
    path = Path(db_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    conn = duckdb.connect(str(path))
    init_schema(conn)
    return conn


def init_schema(conn: duckdb.DuckDBPyConnection) -> None:
    """Create and lightly migrate registry tables."""
    had_legacy_matrix_profiles = _table_exists(conn, "matrix_store_profiles")
    conn.execute(SCHEMA_SQL)
    conn.execute("ALTER TABLE sra_runs ADD COLUMN IF NOT EXISTS status VARCHAR DEFAULT 'added'")
    conn.execute("ALTER TABLE profiles ADD COLUMN IF NOT EXISTS gene_stats_file VARCHAR")
    conn.execute("ALTER TABLE profile_positions ADD COLUMN IF NOT EXISTS ref_base_bitmask UTINYINT")
    conn.execute("ALTER TABLE genome_stats ADD COLUMN IF NOT EXISTS ref_ani DOUBLE")
    conn.execute("ALTER TABLE gene_stats ADD COLUMN IF NOT EXISTS ref_ani DOUBLE")
    conn.execute("ALTER TABLE matrix_stores ADD COLUMN IF NOT EXISTS storage_layout VARCHAR DEFAULT 'dense'")
    conn.execute("ALTER TABLE matrix_stores ADD COLUMN IF NOT EXISTS min_coverage DOUBLE")
    conn.execute("ALTER TABLE matrix_stores ADD COLUMN IF NOT EXISTS min_breadth DOUBLE")
    conn.execute("ALTER TABLE matrix_stores ADD COLUMN IF NOT EXISTS min_ber DOUBLE")
    conn.execute("ALTER TABLE matrix_stores ADD COLUMN IF NOT EXISTS min_sylph_abundance DOUBLE")
    _migrate_profile_counts_to_uint16(conn)
    _normalize_existing_sylph_genomes(conn)
    if had_legacy_matrix_profiles:
        conn.execute(
            """
            INSERT OR IGNORE INTO matrix_store_samples
            SELECT matrix_id, run_id AS sample_id, added_at
            FROM matrix_store_profiles
            """
        )


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
) -> None:
    """Import profile, stat, and abundance files into project tables."""
    if conn.execute("SELECT 1 FROM sra_runs WHERE run_id = ?", [bundle.run_id]).fetchone() is None:
        if not add_run_if_missing:
            raise ValueError(f"Cannot import profile for unknown run_id: {bundle.run_id}")
        add_runs(conn, [bundle.run_id])
    _require_existing_file(bundle.profile_file, "profile_file")
    _require_existing_file(bundle.genome_stats_file, "genome_stats_file")
    _require_existing_file(bundle.sylph_abundance_file, "sylph_abundance_file")
    if bundle.gene_stats_file is not None:
        _require_existing_file(bundle.gene_stats_file, "gene_stats_file")

    sample_id = bundle.run_id
    now = time.time()
    conn.execute("BEGIN TRANSACTION")
    try:
        conn.execute("DELETE FROM profile_positions WHERE sample_id = ?", [sample_id])
        conn.execute("DELETE FROM genome_stats WHERE sample_id = ?", [sample_id])
        conn.execute("DELETE FROM gene_stats WHERE sample_id = ?", [sample_id])
        conn.execute("DELETE FROM sylph_abundance WHERE sample_id = ?", [sample_id])

        _insert_profile_positions(conn, sample_id=sample_id, profile_file=bundle.profile_file)
        _insert_genome_stats(conn, sample_id=sample_id, stats_file=bundle.genome_stats_file)
        if bundle.gene_stats_file is not None:
            _insert_gene_stats(conn, sample_id=sample_id, stats_file=bundle.gene_stats_file)
        _insert_sylph_abundance(conn, sample_id=sample_id, abundance_file=bundle.sylph_abundance_file)

        existing = conn.execute("SELECT created_at FROM samples WHERE sample_id = ?", [sample_id]).fetchone()
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
              (run_id, profile_file, genome_stats_file, gene_stats_file, sylph_abundance_file, created_at, updated_at)
            VALUES (?, ?, ?, ?, ?, ?, ?)
            """,
            [
                bundle.run_id,
                str(bundle.profile_file),
                str(bundle.genome_stats_file),
                str(bundle.gene_stats_file) if bundle.gene_stats_file is not None else None,
                str(bundle.sylph_abundance_file),
                profile_created_at,
                now,
            ],
        )
        conn.execute(
            "UPDATE sra_runs SET status = 'complete', updated_at = ? WHERE run_id = ?",
            [now, bundle.run_id],
        )
        conn.execute("COMMIT")
    except Exception:
        conn.execute("ROLLBACK")
        raise


def add_profiles(
    conn: duckdb.DuckDBPyConnection,
    bundles: list[ProfileBundle],
    *,
    add_runs_if_missing: bool = False,
) -> int:
    """Import completed profile bundles by SRA run ID."""
    for bundle in bundles:
        import_profile_bundle(conn, bundle, add_run_if_missing=add_runs_if_missing)
    return len(bundles)


def list_profiles(conn: duckdb.DuckDBPyConnection) -> list[dict[str, object]]:
    """List completed samples and their source files."""
    return _rows_as_dicts(
        conn.execute(
            """
            SELECT run_id, profile_file, genome_stats_file, gene_stats_file, sylph_abundance_file, created_at, updated_at
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
    progress_callback: ExportProgressCallback | None = None,
) -> list[Path]:
    """Export selected samples from DuckDB into temporary ZipStrain profile parquets."""
    output_dir.mkdir(parents=True, exist_ok=True)
    paths: list[Path] = []
    total = len(sample_ids)
    if progress_callback is not None:
        progress_callback({"phase": "start", "completed": 0, "total": total, "genome": genome or "all"})
    for index, sample_id in enumerate(sample_ids, start=1):
        output_file = output_dir / f"{sample_id}.parquet"
        conditions = ["sample_id = ?"]
        params: list[object] = [sample_id]
        if genome is not None and genome != "all":
            conditions.append("genome = ?")
            params.append(genome)
        where_sql = " AND ".join(conditions)
        arrow_table = conn.execute(
            f"""
            SELECT chrom, genome, pos, COALESCE(gene, 'NA') AS gene, A, T, C, G
            FROM profile_positions
            WHERE {where_sql}
            ORDER BY chrom, pos
            """,
            params,
        ).fetch_arrow_table()
        pl.from_arrow(arrow_table).write_parquet(output_file)
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
    """Record samples materialized into a matrix store."""
    now = time.time()
    inserted = 0
    for sample_id in _normalize_ids(sample_ids):
        row = conn.execute(
            "INSERT OR IGNORE INTO matrix_store_samples VALUES (?, ?, ?) RETURNING sample_id",
            [matrix_id, sample_id, now],
        ).fetchone()
        if row is not None:
            inserted += 1
    return inserted


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


def registry_status(conn: duckdb.DuckDBPyConnection) -> dict[str, int]:
    """Return high-level registry counts."""
    active_runs = conn.execute("SELECT count(*) FROM sra_runs WHERE deleted_at IS NULL").fetchone()[0]
    deleted_runs = conn.execute("SELECT count(*) FROM sra_runs WHERE deleted_at IS NOT NULL").fetchone()[0]
    samples = conn.execute("SELECT count(*) FROM samples WHERE status = 'complete'").fetchone()[0]
    profile_rows = conn.execute("SELECT count(*) FROM profile_positions").fetchone()[0]
    matrices = conn.execute("SELECT count(*) FROM matrix_stores").fetchone()[0]
    compares = conn.execute("SELECT count(*) FROM matrix_compares").fetchone()[0]
    remaining = len(remaining_runs(conn))
    return {
        "active_runs": int(active_runs),
        "deleted_runs": int(deleted_runs),
        "complete_samples": int(samples),
        "profile_rows": int(profile_rows),
        "remaining_profiles": remaining,
        "matrix_stores": int(matrices),
        "matrix_compares": int(compares),
    }


def _insert_profile_positions(conn: duckdb.DuckDBPyConnection, *, sample_id: str, profile_file: Path) -> None:
    df = pl.read_parquet(profile_file)
    required = {"chrom", "pos", "genome", "A", "C", "G", "T"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"profile_file missing required columns: {', '.join(sorted(missing))}")
    if "gene" not in df.columns:
        df = df.with_columns(pl.lit("NA").alias("gene"))
    df = df.select(
        pl.lit(sample_id).alias("sample_id"),
        pl.col("chrom").cast(pl.Utf8),
        pl.col("pos").cast(pl.Int64),
        pl.col("genome").cast(pl.Utf8),
        pl.col("gene").cast(pl.Utf8),
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
               (sample_id, chrom, pos, genome, gene, A, C, G, T, ref_base_bitmask)
               SELECT sample_id, chrom, pos, genome, gene, A, C, G, T, ref_base_bitmask
               FROM _metatrawl_profile_positions"""
        )
    finally:
        conn.unregister("_metatrawl_profile_positions")


def _insert_genome_stats(conn: duckdb.DuckDBPyConnection, *, sample_id: str, stats_file: Path) -> None:
    df = _read_table(stats_file)
    genome_col = _first_existing(df, ["genome", "genome_name", "reference", "accession"])
    if genome_col is None:
        raise ValueError("genome_stats_file missing a genome column")
    df = df.select(
        pl.lit(sample_id).alias("sample_id"),
        pl.col(genome_col).cast(pl.Utf8).alias("genome"),
        _optional_float_expr(df, ["coverage", "cov", "mean_coverage"]).alias("coverage"),
        _optional_float_expr(df, ["breadth", "breadth_coverage", "breadth_cov"]).alias("breadth"),
        _optional_float_expr(df, ["ber", "BER"]).alias("ber"),
        _optional_float_expr(df, ["ref_ani", "reference_ani"]).alias("ref_ani"),
    )
    conn.register("_metatrawl_genome_stats", df)
    try:
        conn.execute(
            """INSERT INTO genome_stats (sample_id, genome, coverage, breadth, ber, ref_ani)
               SELECT sample_id, genome, coverage, breadth, ber, ref_ani FROM _metatrawl_genome_stats"""
        )
    finally:
        conn.unregister("_metatrawl_genome_stats")


def _insert_gene_stats(conn: duckdb.DuckDBPyConnection, *, sample_id: str, stats_file: Path) -> None:
    df = _read_table(stats_file)
    gene_col = _first_existing(df, ["gene", "gene_id"])
    if gene_col is None:
        raise ValueError("gene_stats_file missing a gene column")
    genome_col = _first_existing(df, ["genome", "genome_name", "reference", "accession"])
    df = df.select(
        pl.lit(sample_id).alias("sample_id"),
        (pl.col(genome_col).cast(pl.Utf8) if genome_col else pl.lit(None, dtype=pl.Utf8)).alias("genome"),
        pl.col(gene_col).cast(pl.Utf8).alias("gene"),
        _optional_float_expr(df, ["coverage", "cov", "mean_coverage"]).alias("coverage"),
        _optional_float_expr(df, ["breadth", "breadth_coverage", "breadth_cov"]).alias("breadth"),
        _optional_float_expr(df, ["ber", "BER"]).alias("ber"),
        _optional_float_expr(df, ["ref_ani", "reference_ani"]).alias("ref_ani"),
    )
    conn.register("_metatrawl_gene_stats", df)
    try:
        conn.execute(
            """INSERT INTO gene_stats (sample_id, genome, gene, coverage, breadth, ber, ref_ani)
               SELECT sample_id, genome, gene, coverage, breadth, ber, ref_ani FROM _metatrawl_gene_stats"""
        )
    finally:
        conn.unregister("_metatrawl_gene_stats")


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


def _first_existing(df: pl.DataFrame, names: list[str]) -> str | None:
    columns = set(df.columns)
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


def _normalize_existing_sylph_genomes(conn: duckdb.DuckDBPyConnection) -> None:
    """Replace legacy Sylph database paths with canonical assembly accessions."""
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

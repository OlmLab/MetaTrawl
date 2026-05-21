"""DuckDB registry for MetaTrawl projects."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import time

import duckdb


SCHEMA_SQL = """
CREATE TABLE IF NOT EXISTS sra_runs (
    run_id VARCHAR PRIMARY KEY,
    added_at DOUBLE NOT NULL,
    updated_at DOUBLE NOT NULL,
    deleted_at DOUBLE
);

CREATE TABLE IF NOT EXISTS profiles (
    run_id VARCHAR PRIMARY KEY,
    profile_file VARCHAR NOT NULL,
    genome_stats_file VARCHAR NOT NULL,
    sylph_abundance_file VARCHAR NOT NULL,
    created_at DOUBLE NOT NULL,
    updated_at DOUBLE NOT NULL
);

CREATE TABLE IF NOT EXISTS matrix_stores (
    matrix_id VARCHAR PRIMARY KEY,
    genome VARCHAR NOT NULL,
    matrix_file VARCHAR NOT NULL,
    profile_count BIGINT NOT NULL,
    created_at DOUBLE NOT NULL,
    updated_at DOUBLE NOT NULL
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


@dataclass(frozen=True)
class MatrixStore:
    """Registered matrix store metadata."""

    matrix_id: str
    genome: str
    matrix_file: Path
    profile_count: int


def connect(db_path: str | Path) -> duckdb.DuckDBPyConnection:
    """Open a DuckDB connection and ensure the MetaTrawl schema exists."""
    path = Path(db_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    conn = duckdb.connect(str(path))
    init_schema(conn)
    return conn


def init_schema(conn: duckdb.DuckDBPyConnection) -> None:
    """Create registry tables if they do not already exist."""
    conn.execute(SCHEMA_SQL)


def add_runs(conn: duckdb.DuckDBPyConnection, run_ids: list[str]) -> tuple[int, int]:
    """Add or reactivate SRA run IDs. Returns ``(added, reactivated)``."""
    now = time.time()
    added = 0
    reactivated = 0
    for run_id in _normalize_run_ids(run_ids):
        existing = conn.execute(
            "SELECT deleted_at FROM sra_runs WHERE run_id = ?",
            [run_id],
        ).fetchone()
        if existing is None:
            conn.execute(
                "INSERT INTO sra_runs VALUES (?, ?, ?, NULL)",
                [run_id, now, now],
            )
            added += 1
        elif existing[0] is not None:
            conn.execute(
                "UPDATE sra_runs SET deleted_at = NULL, updated_at = ? WHERE run_id = ?",
                [now, run_id],
            )
            reactivated += 1
    return added, reactivated


def delete_runs(conn: duckdb.DuckDBPyConnection, run_ids: list[str]) -> int:
    """Soft-delete active SRA runs and return the number changed."""
    normalized = _normalize_run_ids(run_ids)
    if not normalized:
        return 0
    now = time.time()
    changed = 0
    for run_id in normalized:
        result = conn.execute(
            """
            UPDATE sra_runs
            SET deleted_at = ?, updated_at = ?
            WHERE run_id = ? AND deleted_at IS NULL
            RETURNING run_id
            """,
            [now, now, run_id],
        ).fetchone()
        if result is not None:
            changed += 1
    return changed


def list_runs(conn: duckdb.DuckDBPyConnection, include_deleted: bool = False) -> list[dict[str, object]]:
    """List registered SRA runs."""
    where = "" if include_deleted else "WHERE deleted_at IS NULL"
    return _rows_as_dicts(
        conn.execute(
            f"""
            SELECT run_id, added_at, updated_at, deleted_at
            FROM sra_runs
            {where}
            ORDER BY run_id
            """
        )
    )


def remaining_runs(conn: duckdb.DuckDBPyConnection) -> list[str]:
    """Return active SRA runs without a complete imported profile bundle."""
    rows = conn.execute(
        """
        SELECT r.run_id
        FROM sra_runs r
        LEFT JOIN profiles p USING (run_id)
        WHERE r.deleted_at IS NULL
          AND (
            p.run_id IS NULL
            OR p.profile_file IS NULL
            OR p.genome_stats_file IS NULL
            OR p.sylph_abundance_file IS NULL
          )
        ORDER BY r.run_id
        """
    ).fetchall()
    return [str(row[0]) for row in rows]


def add_profiles(
    conn: duckdb.DuckDBPyConnection,
    bundles: list[ProfileBundle],
    *,
    add_runs_if_missing: bool = False,
) -> int:
    """Import completed profile bundles by SRA run ID."""
    now = time.time()
    imported = 0
    for bundle in bundles:
        run_exists = conn.execute(
            "SELECT 1 FROM sra_runs WHERE run_id = ?",
            [bundle.run_id],
        ).fetchone()
        if run_exists is None:
            if not add_runs_if_missing:
                raise ValueError(f"Cannot add profile for unknown run_id: {bundle.run_id}")
            add_runs(conn, [bundle.run_id])

        _require_existing_file(bundle.profile_file, "profile_file")
        _require_existing_file(bundle.genome_stats_file, "genome_stats_file")
        _require_existing_file(bundle.sylph_abundance_file, "sylph_abundance_file")

        existing = conn.execute(
            "SELECT created_at FROM profiles WHERE run_id = ?",
            [bundle.run_id],
        ).fetchone()
        created_at = float(existing[0]) if existing is not None else now
        conn.execute(
            """
            INSERT OR REPLACE INTO profiles
            VALUES (?, ?, ?, ?, ?, ?)
            """,
            [
                bundle.run_id,
                str(bundle.profile_file),
                str(bundle.genome_stats_file),
                str(bundle.sylph_abundance_file),
                created_at,
                now,
            ],
        )
        imported += 1
    return imported


def list_profiles(conn: duckdb.DuckDBPyConnection) -> list[dict[str, object]]:
    """List imported profile bundles."""
    return _rows_as_dicts(
        conn.execute(
            """
            SELECT run_id, profile_file, genome_stats_file, sylph_abundance_file, created_at, updated_at
            FROM profiles
            ORDER BY run_id
            """
        )
    )


def complete_profile_paths(conn: duckdb.DuckDBPyConnection) -> list[Path]:
    """Return profile parquet paths for active runs with complete profile bundles."""
    rows = conn.execute(
        """
        SELECT p.profile_file
        FROM profiles p
        JOIN sra_runs r USING (run_id)
        WHERE r.deleted_at IS NULL
        ORDER BY p.run_id
        """
    ).fetchall()
    return [Path(row[0]) for row in rows]


def register_matrix_store(
    conn: duckdb.DuckDBPyConnection,
    *,
    matrix_id: str,
    genome: str,
    matrix_file: Path,
    profile_count: int,
    overwrite: bool = False,
) -> MatrixStore:
    """Record a matrix store in the registry."""
    existing = get_matrix_store(conn, matrix_id)
    if existing is not None and not overwrite:
        raise ValueError(f"Matrix ID already exists: {matrix_id}")
    now = time.time()
    created_at = now
    if existing is not None:
        created_row = conn.execute(
            "SELECT created_at FROM matrix_stores WHERE matrix_id = ?",
            [matrix_id],
        ).fetchone()
        created_at = float(created_row[0])
    conn.execute(
        "INSERT OR REPLACE INTO matrix_stores VALUES (?, ?, ?, ?, ?, ?)",
        [matrix_id, genome, str(matrix_file), profile_count, created_at, now],
    )
    return MatrixStore(matrix_id=matrix_id, genome=genome, matrix_file=matrix_file, profile_count=profile_count)


def get_matrix_store(conn: duckdb.DuckDBPyConnection, matrix_id: str) -> MatrixStore | None:
    """Return a registered matrix store by ID, if present."""
    row = conn.execute(
        "SELECT matrix_id, genome, matrix_file, profile_count FROM matrix_stores WHERE matrix_id = ?",
        [matrix_id],
    ).fetchone()
    if row is None:
        return None
    return MatrixStore(matrix_id=row[0], genome=row[1], matrix_file=Path(row[2]), profile_count=int(row[3]))


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
    existing = conn.execute(
        "SELECT created_at FROM matrix_compares WHERE compare_id = ?",
        [compare_id],
    ).fetchone()
    created_at = float(existing[0]) if existing is not None else now
    conn.execute(
        "INSERT OR REPLACE INTO matrix_compares VALUES (?, ?, ?, ?, ?, ?)",
        [compare_id, matrix_id, str(compare_db_file), calculate, created_at, now],
    )


def registry_status(conn: duckdb.DuckDBPyConnection) -> dict[str, int]:
    """Return high-level registry counts."""
    active_runs = conn.execute("SELECT count(*) FROM sra_runs WHERE deleted_at IS NULL").fetchone()[0]
    deleted_runs = conn.execute("SELECT count(*) FROM sra_runs WHERE deleted_at IS NOT NULL").fetchone()[0]
    profiles = conn.execute("SELECT count(*) FROM profiles").fetchone()[0]
    matrices = conn.execute("SELECT count(*) FROM matrix_stores").fetchone()[0]
    compares = conn.execute("SELECT count(*) FROM matrix_compares").fetchone()[0]
    remaining = len(remaining_runs(conn))
    return {
        "active_runs": int(active_runs),
        "deleted_runs": int(deleted_runs),
        "profiles": int(profiles),
        "remaining_profiles": remaining,
        "matrix_stores": int(matrices),
        "matrix_compares": int(compares),
    }


def _normalize_run_ids(run_ids: list[str]) -> list[str]:
    seen: set[str] = set()
    normalized: list[str] = []
    for run_id in run_ids:
        clean = str(run_id).strip()
        if clean and clean not in seen:
            normalized.append(clean)
            seen.add(clean)
    return normalized


def _require_existing_file(path: Path, column: str) -> None:
    if not path.exists():
        raise FileNotFoundError(f"{column} does not exist: {path}")


def _rows_as_dicts(result) -> list[dict[str, object]]:
    columns = [desc[0] for desc in result.description]
    return [dict(zip(columns, row)) for row in result.fetchall()]

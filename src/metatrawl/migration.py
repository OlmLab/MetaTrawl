"""Resumable migration from full profile rows to allele-mask storage."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import time
from typing import Callable

import duckdb

from metatrawl import allele_mask
from metatrawl import db


MigrationProgressCallback = Callable[[dict[str, object]], None]
_COPY_BATCH_ROWS = 262_144
_REGISTRY_TABLES = (
    "sra_runs",
    "samples",
    "profiles",
    "genome_stats",
    "gene_stats",
    "sylph_abundance",
    "cache_genomes",
    "matrix_stores",
    "matrix_store_samples",
    "matrix_compares",
)


@dataclass(frozen=True)
class AlleleMaskMigrationSummary:
    """Result of one resumable database migration invocation."""

    total_samples: int
    completed_samples: int
    migrated_samples: int
    failed_samples: int


def migrate_full_database(
    *,
    source_db: Path,
    output_db: Path,
    cache_dir: Path,
    min_cov: int,
    progress_callback: MigrationProgressCallback | None = None,
) -> AlleleMaskMigrationSummary:
    """Create or resume an allele-mask database without modifying the source."""
    source_db = Path(source_db).expanduser().resolve()
    output_db = Path(output_db).expanduser().resolve()
    cache_dir = Path(cache_dir).expanduser().resolve()
    if source_db == output_db:
        raise ValueError("Source and output databases must be different files.")
    if not source_db.is_file():
        raise FileNotFoundError(f"Source database does not exist: {source_db}")
    if not cache_dir.is_dir():
        raise FileNotFoundError(f"Genome cache directory does not exist: {cache_dir}")
    if min_cov < 1:
        raise ValueError("min_cov must be >= 1.")

    output_existed = output_db.exists()
    output_db.parent.mkdir(parents=True, exist_ok=True)
    source = duckdb.connect(str(source_db), read_only=True)
    target = db.connect(output_db)
    try:
        source_storage = db.profile_storage_config(source)
        if source_storage.mode != allele_mask.PROFILE_STORAGE_FULL:
            raise ValueError(
                "Migration source must use full profile storage; "
                f"found {source_storage.mode}."
            )
        if output_existed:
            target_storage = db.profile_storage_config(target)
            state_rows = int(
                target.execute(
                    "SELECT count(*) FROM allele_mask_migration_state"
                ).fetchone()[0]
            )
            if _is_empty_uninitialized_target(target):
                _initialize_target(
                    source,
                    target,
                    min_cov=min_cov,
                    reset=True,
                    progress_callback=progress_callback,
                )
            elif (
                target_storage.mode != allele_mask.PROFILE_STORAGE_ALLELE_MASK
                or target_storage.min_cov != min_cov
            ):
                raise ValueError(
                    "Existing migration output has a different profile-storage "
                    "contract. Use the same --min-cov or a new output file."
                )
            elif state_rows == 0 and int(
                source.execute(
                    "SELECT count(*) FROM samples WHERE status = 'complete'"
                ).fetchone()[0]
            ):
                raise ValueError(
                    "Output database uses allele-mask storage but has no migration "
                    "checkpoint. Refusing to overwrite it; choose a new output file."
                )
        else:
            _initialize_target(
                source,
                target,
                min_cov=min_cov,
                reset=False,
                progress_callback=progress_callback,
            )

        pending = [
            str(row[0])
            for row in target.execute(
                """
                SELECT sample_id
                FROM allele_mask_migration_state
                WHERE status != 'done'
                ORDER BY sample_id
                """
            ).fetchall()
        ]
        total = int(
            target.execute(
                "SELECT count(*) FROM allele_mask_migration_state"
            ).fetchone()[0]
        )
        completed_before = total - len(pending)
        _emit(
            progress_callback,
            phase="start",
            total=total,
            completed=completed_before,
            remaining=len(pending),
        )
        migrated = 0
        failed = 0
        for index, sample_id in enumerate(pending, start=completed_before + 1):
            _emit(
                progress_callback,
                phase="sample-start",
                sample_id=sample_id,
                completed=index - 1,
                total=total,
            )
            try:
                reader = source.execute(
                    """
                    SELECT chrom, genome, pos, A, C, G, T
                    FROM profile_positions
                    WHERE sample_id = ?
                    """,
                    [sample_id],
                ).fetch_record_batch(allele_mask.PROFILE_BATCH_ROWS)
                target.execute("BEGIN TRANSACTION")
                target.execute(
                    "DELETE FROM allele_mask_profile_blocks WHERE sample_id = ?",
                    [sample_id],
                )
                summary = allele_mask.store_profile_batches(
                    target,
                    sample_id=sample_id,
                    batches=reader,
                    min_cov=min_cov,
                    cache_dir=cache_dir,
                )
                target.execute(
                    """
                    UPDATE allele_mask_migration_state
                    SET status = 'done', source_rows = ?, compressed_blocks = ?,
                        covered_positions = ?, error = NULL, updated_at = ?
                    WHERE sample_id = ?
                    """,
                    [
                        summary.source_rows,
                        summary.blocks,
                        summary.covered_positions,
                        time.time(),
                        sample_id,
                    ],
                )
                target.execute(
                    "UPDATE samples SET status = 'complete', updated_at = ? WHERE sample_id = ?",
                    [time.time(), sample_id],
                )
                target.execute("COMMIT")
            except Exception as exc:
                try:
                    target.execute("ROLLBACK")
                except duckdb.TransactionException:
                    pass
                target.execute(
                    """
                    UPDATE allele_mask_migration_state
                    SET status = 'failed', error = ?, updated_at = ?
                    WHERE sample_id = ?
                    """,
                    [str(exc), time.time(), sample_id],
                )
                failed += 1
                _emit(
                    progress_callback,
                    phase="sample-failed",
                    sample_id=sample_id,
                    completed=index - 1,
                    total=total,
                    error=str(exc),
                )
                continue
            migrated += 1
            _emit(
                progress_callback,
                phase="sample-done",
                sample_id=sample_id,
                completed=index,
                total=total,
                source_rows=summary.source_rows,
                blocks=summary.blocks,
                covered_positions=summary.covered_positions,
            )

        completed = int(
            target.execute(
                "SELECT count(*) FROM allele_mask_migration_state WHERE status = 'done'"
            ).fetchone()[0]
        )
        _emit(
            progress_callback,
            phase="done",
            total=total,
            completed=completed,
            migrated=migrated,
            failed=failed,
        )
        return AlleleMaskMigrationSummary(
            total_samples=total,
            completed_samples=completed,
            migrated_samples=migrated,
            failed_samples=failed,
        )
    finally:
        target.close()
        source.close()


def _initialize_target(
    source: duckdb.DuckDBPyConnection,
    target: duckdb.DuckDBPyConnection,
    *,
    min_cov: int,
    reset: bool,
    progress_callback: MigrationProgressCallback | None,
) -> None:
    """Atomically publish copied metadata and initial per-sample checkpoints."""
    target.execute("BEGIN TRANSACTION")
    try:
        if reset:
            target.execute("DELETE FROM allele_mask_profile_blocks")
            target.execute("DELETE FROM allele_mask_reference_segments")
            target.execute("DELETE FROM allele_mask_migration_state")
            for table in reversed(_REGISTRY_TABLES):
                target.execute(f"DELETE FROM {_quote_identifier(table)}")
        db.configure_profile_storage(
            target,
            mode=allele_mask.PROFILE_STORAGE_ALLELE_MASK,
            min_cov=min_cov,
        )
        _copy_registry_tables(
            source,
            target,
            progress_callback=progress_callback,
        )
        target.execute(
            """
            UPDATE profiles
            SET profile_storage_mode = 'allele-mask', profile_min_cov = ?
            """,
            [min_cov],
        )
        target.execute(
            """
            INSERT INTO allele_mask_migration_state
              (sample_id, status, source_rows, compressed_blocks,
               covered_positions, error, updated_at)
            SELECT sample_id, 'pending', NULL, NULL, NULL, NULL, ?
            FROM samples
            WHERE status = 'complete'
            """,
            [time.time()],
        )
        target.execute(
            """
            UPDATE samples
            SET status = 'migrating'
            WHERE sample_id IN (
                SELECT sample_id FROM allele_mask_migration_state
            )
            """
        )
        target.execute("COMMIT")
    except Exception:
        target.execute("ROLLBACK")
        raise


def _is_empty_uninitialized_target(
    target: duckdb.DuckDBPyConnection,
) -> bool:
    storage = db.profile_storage_config(target)
    if storage.mode != allele_mask.PROFILE_STORAGE_FULL:
        return False
    table_counts = [
        int(target.execute(f"SELECT count(*) FROM {_quote_identifier(table)}").fetchone()[0])
        for table in _REGISTRY_TABLES
    ]
    compact_count = int(
        target.execute(
            """
            SELECT
              (SELECT count(*) FROM allele_mask_reference_segments)
              + (SELECT count(*) FROM allele_mask_profile_blocks)
              + (SELECT count(*) FROM allele_mask_migration_state)
            """
        ).fetchone()[0]
    )
    return not any(table_counts) and compact_count == 0


def _copy_registry_tables(
    source: duckdb.DuckDBPyConnection,
    target: duckdb.DuckDBPyConnection,
    *,
    progress_callback: MigrationProgressCallback | None,
) -> None:
    for table in _REGISTRY_TABLES:
        if not _table_exists(source, table):
            continue
        source_columns = _table_columns(source, table)
        target_columns = set(_table_columns(target, table))
        columns = [column for column in source_columns if column in target_columns]
        if not columns:
            continue
        _emit(progress_callback, phase="copy-table", table=table)
        quoted = ", ".join(_quote_identifier(column) for column in columns)
        reader = source.execute(
            f"SELECT {quoted} FROM {_quote_identifier(table)}"
        ).fetch_record_batch(_COPY_BATCH_ROWS)
        for batch in reader:
            target.register("_metatrawl_migration_batch", batch)
            try:
                target.execute(
                    f"""
                    INSERT INTO {_quote_identifier(table)} ({quoted})
                    SELECT {quoted}
                    FROM _metatrawl_migration_batch
                    """
                )
            finally:
                target.unregister("_metatrawl_migration_batch")


def _table_exists(conn: duckdb.DuckDBPyConnection, table: str) -> bool:
    return (
        conn.execute(
            """
            SELECT 1
            FROM information_schema.tables
            WHERE table_schema = 'main' AND table_name = ?
            """,
            [table],
        ).fetchone()
        is not None
    )


def _table_columns(conn: duckdb.DuckDBPyConnection, table: str) -> list[str]:
    return [
        str(row[1])
        for row in conn.execute(
            f"PRAGMA table_info({_quote_string(table)})"
        ).fetchall()
    ]


def _quote_identifier(value: str) -> str:
    return '"' + value.replace('"', '""') + '"'


def _quote_string(value: str) -> str:
    return "'" + value.replace("'", "''") + "'"


def _emit(
    callback: MigrationProgressCallback | None,
    **event: object,
) -> None:
    if callback is not None:
        callback(event)

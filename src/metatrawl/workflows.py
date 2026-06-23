"""Workflow adapters for MetaTrawl."""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
import re
import shlex
import shutil
import subprocess
from tempfile import TemporaryDirectory
from typing import Callable

import polars as pl

from metatrawl import cache
from metatrawl import db
from metatrawl import healthcheck
from metatrawl.config import MatrixCompareConfig, WorkflowConfig
from metatrawl.execution import WorkflowRuntime
from metatrawl.logging import ThrottledMatrixLogger, WorkflowLogger


ACCESSION_PATTERN = re.compile(r"\b((?:GC[AF])_\d+(?:\.\d+)?)", re.IGNORECASE)
SYLPH_GENOME_COLUMNS = ["Genome_file", "genome_file", "genome", "Genome", "Reference", "reference", "target"]
SYLPH_ABUNDANCE_COLUMNS = ["eff_cov", "genome_cov", "abundance", "relative_abundance", "rel_abundance", "Taxonomic_abundance"]


@dataclass(frozen=True)
class SyncSummary:
    """High-level result for one MetaTrawl sync run."""

    requested: int
    imported: int
    skipped: int
    failed: int
    cleaned_files: int


@dataclass(frozen=True)
class MatrixFileContext:
    """Metadata read directly from a ZipStrain HDF5 matrix store."""

    matrix_file: Path
    genome: str
    storage_layout: str
    sample_ids: tuple[str, ...]
    filters: db.MatrixFilters


@dataclass(frozen=True)
class MatrixSyncBuildSummary:
    """Result from converging all per-genome matrix files."""

    genomes: int
    built: int
    appended: int
    up_to_date: int
    skipped: int
    failed: int


@dataclass(frozen=True)
class MatrixSyncCompareSummary:
    """Result from converging compare databases for matrix files."""

    matrices: int
    compared: int
    failed: int


def build_matrix_from_database(
    conn,
    *,
    output_file: Path,
    genome: str,
    bed_file: Path,
    stb_file: Path,
    filters: db.MatrixFilters,
    overwrite: bool = False,
    matrix_id: str | None = None,
    count_dtype: str = "uint16",
    gene_range_table: Path | None = None,
    memory_limit_gb: float = 16.0,
    export_batch_mb: float = 128.0,
    sparse: bool = False,
    logger: WorkflowLogger | None = None,
) -> db.MatrixStore:
    """Build a ZipStrain matrix store from DuckDB profile rows."""
    logger = logger or WorkflowLogger()
    sample_ids = db.eligible_sample_ids(conn, genome=genome, filters=filters)
    if not sample_ids:
        raise ValueError("No complete samples passed the matrix build filters.")

    output_file = Path(output_file)
    if output_file.exists() and not overwrite:
        raise FileExistsError(f"Matrix output already exists: {output_file}")
    output_file.parent.mkdir(parents=True, exist_ok=True)

    chosen_matrix_id = matrix_id or output_file.stem

    from zipstrain import matrix_pairs as mp

    logger.emit(step="matrix-build", status="exporting-profiles", samples=len(sample_ids), genome=genome)
    with TemporaryDirectory(prefix="metatrawl_matrix_profiles_") as tmp_dir:
        profile_dir = Path(tmp_dir)
        db.export_profile_parquets(
            conn,
            sample_ids=sample_ids,
            output_dir=profile_dir,
            genome=None if genome == "all" else genome,
            progress_callback=ThrottledMatrixLogger("METATRAWL-EXPORT", stored_rows=False),
        )
        logger.emit(step="matrix-build", status="exported-profiles", samples=len(sample_ids), genome=genome)
        logger.emit(step="matrix-build", status="building", samples=len(sample_ids), sparse=sparse)
        mp.build_matrix_hdf5(
            profile_dir=profile_dir,
            output_file=output_file,
            genome=genome,
            bed_file=bed_file,
            stb_file=stb_file,
            gene_range_table=gene_range_table,
            count_dtype=count_dtype,
            memory_limit_gb=memory_limit_gb,
            export_batch_mb=export_batch_mb,
            sparse=sparse,
            progress_callback=ThrottledMatrixLogger("MATRIX-BUILD"),
        )

    store = db.register_matrix_store(
        conn,
        matrix_id=chosen_matrix_id,
        genome=genome,
        matrix_file=output_file,
        profile_count=len(sample_ids),
        storage_layout="sparse" if sparse else "dense",
        sample_ids=sample_ids,
        filters=filters,
        overwrite=True,
    )
    _write_matrix_hdf5_metatrawl_metadata(output_file, filters=filters)
    logger.emit(step="matrix-build", status="done", matrix_id=chosen_matrix_id, samples=len(sample_ids))
    return store


def append_matrix_from_database(
    conn,
    *,
    matrix_file: Path,
    matrix_id: str | None = None,
    memory_limit_gb: float = 16.0,
    export_batch_mb: float = 128.0,
    logger: WorkflowLogger | None = None,
) -> int:
    """Append newly imported complete samples into a matrix store HDF5 file."""
    logger = logger or WorkflowLogger()
    context = read_matrix_file_context(matrix_file)
    matrix_label = matrix_id or str(context.matrix_file)
    sample_ids = (
        [
            sample_id
            for sample_id in db.completed_sample_ids(conn)
            if sample_id not in set(context.sample_ids)
        ]
        if context.genome == "all"
        else db.eligible_unmaterialized_sample_ids(
            conn,
            matrix_id=None,
            existing_sample_ids=list(context.sample_ids),
            genome=context.genome,
            filters=context.filters,
        )
    )
    if not sample_ids:
        raise ValueError(f"No new complete samples are available to append to matrix: {context.matrix_file}")

    from zipstrain import matrix_pairs as mp

    logger.emit(step="matrix-append", status="exporting-profiles", matrix=matrix_label, samples=len(sample_ids))
    with TemporaryDirectory(prefix="metatrawl_matrix_append_") as tmp_dir:
        profile_dir = Path(tmp_dir)
        db.export_profile_parquets(
            conn,
            sample_ids=sample_ids,
            output_dir=profile_dir,
            genome=None if context.genome == "all" else context.genome,
            progress_callback=ThrottledMatrixLogger("METATRAWL-EXPORT", stored_rows=False),
        )
        logger.emit(step="matrix-append", status="appending", matrix=matrix_label, samples=len(sample_ids))
        mp.append_matrix_hdf5(
            profile_dir=profile_dir,
            matrix_hdf5_file=context.matrix_file,
            memory_limit_gb=memory_limit_gb,
            export_batch_mb=export_batch_mb,
            progress_callback=ThrottledMatrixLogger("MATRIX-APPEND"),
        )
    if matrix_id is not None and db.get_matrix_store(conn, matrix_id) is not None:
        db.add_matrix_store_samples(conn, matrix_id=matrix_id, sample_ids=sample_ids)
        db.update_matrix_profile_count(conn, matrix_id=matrix_id)
    logger.emit(step="matrix-append", status="done", matrix=matrix_label, samples=len(sample_ids))
    return len(sample_ids)


def sync_build_matrices(
    conn,
    *,
    matrix_dir: Path,
    genomes: list[str] | None = None,
    bed_file: Path | None = None,
    stb_file: Path | None = None,
    bed_dir: Path | None = None,
    stb_dir: Path | None = None,
    gene_range_table: Path | None = None,
    gene_range_dir: Path | None = None,
    filters: db.MatrixFilters | None = None,
    sparse: bool = False,
    memory_limit_gb: float = 16.0,
    export_batch_mb: float = 128.0,
    logger: WorkflowLogger | None = None,
) -> MatrixSyncBuildSummary:
    """Build or append one matrix per genome represented in complete samples."""
    logger = logger or WorkflowLogger()
    filters = filters or db.MatrixFilters()
    matrix_dir = Path(matrix_dir)
    matrix_dir.mkdir(parents=True, exist_ok=True)
    selected_genomes = _dedupe_preserving_order(genomes) if genomes is not None else db.genomes_with_complete_samples(conn)
    logger.emit(step="matrix-sync-build", status="start", genomes=len(selected_genomes), matrix_dir=matrix_dir)

    built = appended = up_to_date = skipped = failed = 0
    for genome in selected_genomes:
        matrix_file = matrix_dir / f"{_safe_file_stem(genome)}.h5"
        logger.emit(step="matrix-sync-build", status="genome-start", genome=genome, matrix=matrix_file)
        try:
            if matrix_file.exists():
                try:
                    appended_count = append_matrix_from_database(
                        conn,
                        matrix_file=matrix_file,
                        memory_limit_gb=memory_limit_gb,
                        export_batch_mb=export_batch_mb,
                        logger=logger,
                    )
                except ValueError as exc:
                    if "No new complete samples" not in str(exc):
                        raise
                    up_to_date += 1
                    logger.emit(step="matrix-sync-build", status="up-to-date", genome=genome, matrix=matrix_file)
                else:
                    appended += 1
                    logger.emit(
                        step="matrix-sync-build",
                        status="appended",
                        genome=genome,
                        matrix=matrix_file,
                        samples=appended_count,
                    )
            else:
                store = build_matrix_from_database(
                    conn,
                    output_file=matrix_file,
                    genome=genome,
                    bed_file=_per_genome_file_or_fallback(
                        genome=genome,
                        directory=bed_dir,
                        suffix=".bed",
                        fallback=bed_file,
                    ),
                    stb_file=_per_genome_file_or_fallback(
                        genome=genome,
                        directory=stb_dir,
                        suffix=".stb",
                        fallback=stb_file,
                    ),
                    gene_range_table=_gene_range_table_for_genome(
                        genome=genome,
                        gene_range_dir=gene_range_dir,
                        fallback=gene_range_table,
                    ),
                    filters=filters,
                    sparse=sparse,
                    memory_limit_gb=memory_limit_gb,
                    export_batch_mb=export_batch_mb,
                    logger=logger,
                )
                built += 1
                logger.emit(
                    step="matrix-sync-build",
                    status="built",
                    genome=genome,
                    matrix=matrix_file,
                    samples=store.profile_count,
                )
        except ValueError as exc:
            if "No complete samples passed" in str(exc):
                skipped += 1
                logger.emit(step="matrix-sync-build", status="skipped", genome=genome, matrix=matrix_file, error=exc)
                continue
            failed += 1
            logger.emit(step="matrix-sync-build", status="failed", genome=genome, matrix=matrix_file, error=exc)
        except Exception as exc:
            failed += 1
            logger.emit(step="matrix-sync-build", status="failed", genome=genome, matrix=matrix_file, error=exc)

    logger.emit(
        step="matrix-sync-build",
        status="done",
        genomes=len(selected_genomes),
        built=built,
        appended=appended,
        up_to_date=up_to_date,
        skipped=skipped,
        failed=failed,
    )
    return MatrixSyncBuildSummary(
        genomes=len(selected_genomes),
        built=built,
        appended=appended,
        up_to_date=up_to_date,
        skipped=skipped,
        failed=failed,
    )


def sync_compare_matrices(
    conn,
    *,
    matrix_dir: Path,
    compare_dir: Path,
    calculate: str = "all",
    genome: str = "all",
    backend: str = "numpy",
    memory_limit_gb: float = 16.0,
    compare_config: MatrixCompareConfig | None = None,
    logger: WorkflowLogger | None = None,
) -> MatrixSyncCompareSummary:
    """Run resumable compare for every HDF5 matrix file in a directory."""
    logger = logger or WorkflowLogger()
    matrix_files = discover_matrix_files(matrix_dir)
    compare_dir = Path(compare_dir)
    compare_dir.mkdir(parents=True, exist_ok=True)
    logger.emit(step="matrix-sync-compare", status="start", matrices=len(matrix_files), matrix_dir=matrix_dir)

    compared = failed = 0
    for matrix_file in matrix_files:
        output_file = compare_dir / f"{matrix_file.stem}.duckdb"
        logger.emit(step="matrix-sync-compare", status="matrix-start", matrix=matrix_file, output=output_file)
        try:
            compare_matrix_file(
                conn,
                matrix_file=matrix_file,
                output_file=output_file,
                calculate=calculate,
                genome=genome,
                backend=backend,
                memory_limit_gb=memory_limit_gb,
                compare_config=compare_config,
                logger=logger,
            )
            compared += 1
            logger.emit(step="matrix-sync-compare", status="compared", matrix=matrix_file, output=output_file)
        except Exception as exc:
            failed += 1
            logger.emit(step="matrix-sync-compare", status="failed", matrix=matrix_file, output=output_file, error=exc)

    logger.emit(step="matrix-sync-compare", status="done", matrices=len(matrix_files), compared=compared, failed=failed)
    return MatrixSyncCompareSummary(matrices=len(matrix_files), compared=compared, failed=failed)


def discover_matrix_files(matrix_dir: Path) -> list[Path]:
    """Return HDF5-like matrix files from a matrix directory."""
    matrix_dir = Path(matrix_dir)
    if not matrix_dir.is_dir():
        raise FileNotFoundError(f"Matrix directory does not exist: {matrix_dir}")
    suffixes = {".h5", ".hdf5", ".hd5"}
    return sorted(path for path in matrix_dir.iterdir() if path.is_file() and path.suffix.lower() in suffixes)


def _safe_file_stem(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", value).strip("_") or "matrix"


def _dedupe_preserving_order(values: list[str]) -> list[str]:
    seen: set[str] = set()
    result: list[str] = []
    for value in values:
        clean = value.strip()
        if clean and clean not in seen:
            seen.add(clean)
            result.append(clean)
    return result


def _per_genome_file_or_fallback(*, genome: str, directory: Path | None, suffix: str, fallback: Path | None) -> Path:
    if directory is None:
        if fallback is None:
            raise FileNotFoundError(f"No per-genome {suffix} directory or fallback file was provided for genome: {genome}")
        return fallback
    candidate = Path(directory) / f"{_safe_file_stem(genome)}{suffix}"
    if candidate.exists():
        return candidate
    if fallback is not None:
        return fallback
    raise FileNotFoundError(f"Missing per-genome matrix requirement file: {candidate}")


def _gene_range_table_for_genome(*, genome: str, gene_range_dir: Path | None, fallback: Path | None) -> Path | None:
    if gene_range_dir is None:
        return fallback
    candidate = Path(gene_range_dir) / f"{_safe_file_stem(genome)}.gene_ranges.tsv"
    return candidate if candidate.exists() else fallback


def resolve_matrix_store(conn, *, matrix_id: str | None = None, matrix_file: Path | None = None) -> db.MatrixStore:
    """Resolve a matrix store from either its registry ID or its file path."""
    if (matrix_id is None) == (matrix_file is None):
        raise ValueError("Provide exactly one of matrix_id or matrix_file.")
    if matrix_id is not None:
        store = db.get_matrix_store(conn, matrix_id)
        if store is None:
            raise ValueError(f"Unknown matrix ID: {matrix_id}")
        return store
    assert matrix_file is not None
    store = db.get_matrix_store_by_file(conn, matrix_file)
    if store is None:
        raise ValueError(
            f"Matrix file is not registered in MetaTrawl: {matrix_file}. "
            "Build it with `metatrawl matrix build` first."
        )
    return store


def resolve_matrix_file(conn, *, matrix_id: str | None = None, matrix_file: Path | None = None) -> tuple[Path, str | None]:
    """Resolve a matrix file from either a path or a legacy registry ID."""
    if (matrix_id is None) == (matrix_file is None):
        raise ValueError("Provide exactly one of matrix_id or matrix_file.")
    if matrix_file is not None:
        return Path(matrix_file), None
    assert matrix_id is not None
    store = db.get_matrix_store(conn, matrix_id)
    if store is None:
        raise ValueError(f"Unknown matrix ID: {matrix_id}")
    return store.matrix_file, store.matrix_id


def compare_matrix_from_registry(
    conn,
    *,
    matrix_id: str,
    output_file: Path,
    calculate: str = "all",
    genome: str = "all",
    backend: str = "numpy",
    memory_limit_gb: float = 16.0,
    logger: WorkflowLogger | None = None,
) -> str:
    """Run ZipStrain matrix compare for a registered matrix store."""
    logger = logger or WorkflowLogger()
    store = db.get_matrix_store(conn, matrix_id)
    if store is None:
        raise ValueError(f"Unknown matrix ID: {matrix_id}")
    output_file = Path(output_file)
    output_file.parent.mkdir(parents=True, exist_ok=True)

    from zipstrain import matrix_pairs as mp

    logger.emit(step="matrix-compare", status="start", matrix_id=matrix_id, calculate=calculate)
    mp.matrix_compare(
        matrix_db_file=store.matrix_file,
        output_file=output_file,
        genome=genome,
        memory_limit_gb=memory_limit_gb,
        backend=backend,
        calculate=calculate,
    )
    compare_id = output_file.stem
    db.register_matrix_compare(
        conn,
        compare_id=compare_id,
        matrix_id=matrix_id,
        compare_db_file=output_file,
        calculate=calculate,
    )
    logger.emit(step="matrix-compare", status="done", matrix_id=matrix_id, compare_id=compare_id)
    return compare_id


def compare_matrix_file(
    conn,
    *,
    matrix_file: Path,
    output_file: Path,
    calculate: str = "all",
    genome: str = "all",
    backend: str = "numpy",
    memory_limit_gb: float = 16.0,
    matrix_id: str | None = None,
    compare_config: MatrixCompareConfig | None = None,
    logger: WorkflowLogger | None = None,
) -> str:
    """Run ZipStrain matrix compare using the HDF5 file as the durable handle."""
    logger = logger or WorkflowLogger()
    context = read_matrix_file_context(matrix_file)
    output_file = Path(output_file)
    output_file.parent.mkdir(parents=True, exist_ok=True)

    from zipstrain import matrix_pairs as mp

    logger.emit(step="matrix-compare", status="start", matrix=context.matrix_file, calculate=calculate)
    compare_kwargs = compare_config.kwargs() if compare_config is not None else {}
    mp.matrix_compare(
        matrix_db_file=context.matrix_file,
        output_file=output_file,
        genome=genome,
        memory_limit_gb=memory_limit_gb,
        backend=backend,
        calculate=calculate,
        **compare_kwargs,
    )
    compare_id = output_file.stem
    db.register_matrix_compare(
        conn,
        compare_id=compare_id,
        matrix_id=matrix_id or str(context.matrix_file),
        compare_db_file=output_file,
        calculate=calculate,
    )
    logger.emit(step="matrix-compare", status="done", matrix=context.matrix_file, compare_id=compare_id)
    return compare_id


def read_matrix_file_context(matrix_file: Path) -> MatrixFileContext:
    """Read matrix metadata and materialized samples from a ZipStrain HDF5 store."""
    path = Path(matrix_file).expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError(f"Matrix store file does not exist: {path}")
    h5py = _import_h5py()
    with h5py.File(str(path), "r") as handle:
        if "metadata" not in handle:
            raise ValueError(f"Matrix store is missing metadata group: {path}")
        metadata = {str(key): str(value) for key, value in handle["metadata"].attrs.items()}
        genome = metadata.get("genome_scope", "all")
        layout = metadata.get("layout", "dense")
        if "samples" not in handle or "sample_name" not in handle["samples"]:
            sample_ids: tuple[str, ...] = ()
        else:
            sample_ids = tuple(_decode_hdf5_value(value) for value in handle["samples"]["sample_name"][()])
        filters = db.MatrixFilters(
            min_coverage=_optional_float_metadata(metadata, "metatrawl_min_coverage"),
            min_breadth=_optional_float_metadata(metadata, "metatrawl_min_breadth"),
            min_ber=_optional_float_metadata(metadata, "metatrawl_min_ber"),
            min_sylph_abundance=_optional_float_metadata(metadata, "metatrawl_min_sylph_abundance"),
        )
    return MatrixFileContext(
        matrix_file=path,
        genome=genome,
        storage_layout=layout,
        sample_ids=sample_ids,
        filters=filters,
    )


def _write_matrix_hdf5_metatrawl_metadata(matrix_file: Path, *, filters: db.MatrixFilters) -> None:
    h5py = _import_h5py()
    with h5py.File(str(Path(matrix_file).expanduser().resolve()), "r+") as handle:
        metadata = handle.require_group("metadata")
        metadata.attrs["metatrawl_min_coverage"] = _metadata_float(filters.min_coverage)
        metadata.attrs["metatrawl_min_breadth"] = _metadata_float(filters.min_breadth)
        metadata.attrs["metatrawl_min_ber"] = _metadata_float(filters.min_ber)
        metadata.attrs["metatrawl_min_sylph_abundance"] = _metadata_float(filters.min_sylph_abundance)


def _import_h5py():
    try:
        import h5py  # type: ignore
    except ImportError as exc:
        raise RuntimeError("HDF5 matrix operations require h5py. Install MetaTrawl with matrix dependencies.") from exc
    return h5py


def _decode_hdf5_value(value: object) -> str:
    if isinstance(value, bytes):
        return value.decode("utf-8")
    return str(value)


def _metadata_float(value: float | None) -> str:
    return "" if value is None else str(value)


def _optional_float_metadata(metadata: dict[str, str], key: str) -> float | None:
    value = metadata.get(key)
    if value in (None, ""):
        return None
    return float(value)


def profile_sra_runs(
    *,
    run_ids: list[str],
    db_file: Path,
    cache_dir: Path,
    scratch_dir: Path,
    sylph_db: Path | None = None,
    output_dir: Path | None = None,
    accessions_dir: Path | None = None,
    threads: int = 8,
    workflow_config: WorkflowConfig | None = None,
    logger: WorkflowLogger | None = None,
    raise_on_error: bool = True,
    completion_callback: Callable[[str], None] | None = None,
) -> dict[str, str]:
    """Profile SRA runs, retaining failed stage outputs for resumable retries.

    ``completion_callback`` runs serially in the coordinator thread as soon as a
    sample finishes. Sync uses it to import each sample before its scratch data is
    removed, avoiding concurrent DuckDB writers and batch-wide import delays.
    """
    logger = logger or WorkflowLogger()
    workflow_config = workflow_config or WorkflowConfig.legacy(threads=threads, sample_count=len(run_ids))
    runtime = WorkflowRuntime(workflow_config, state_dir=Path(scratch_dir), logger=logger, runner=subprocess.run)
    cache_manager = cache.GenomeCache(
        cache_dir,
        logger=logger,
        downloader=lambda accession, output: cache.download_genome_with_datasets(
            accession,
            output,
            command_runner=lambda command: runtime.run("genome_download", command, sample=accession),
        ),
        prodigal_runner=lambda genome, output: cache.run_prodigal_gene_fasta(
            genome,
            output,
            command_runner=lambda command: runtime.run("prodigal", command, sample=genome.stem),
        ),
        download_workers=workflow_config.stage("genome_download").workers,
        prodigal_workers=workflow_config.stage("prodigal").workers,
    )
    sylph_db = _resolve_existing_sylph_db(sylph_db) if sylph_db is not None else None
    max_workers = max(1, min(len(run_ids), workflow_config.sample_workers))
    logger.emit(step="profile-sra", status="start", samples=len(run_ids), workers=max_workers, db=db_file)
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(
                _profile_one_sra_run,
                run_id=run_id,
                cache_manager=cache_manager,
                scratch_dir=Path(scratch_dir),
                sylph_db=sylph_db,
                output_dir=Path(output_dir) if output_dir is not None else None,
                accessions_dir=Path(accessions_dir) if accessions_dir is not None else None,
                runtime=runtime,
                logger=logger,
            ): run_id
            for run_id in run_ids
        }
        failures: dict[str, str] = {}
        for future in as_completed(futures):
            run_id = futures[future]
            try:
                future.result()
                if completion_callback is not None:
                    completion_callback(run_id)
                sample_scratch = Path(scratch_dir) / run_id
                if sample_scratch.exists():
                    shutil.rmtree(sample_scratch)
                    logger.emit(sample=run_id, step="cleanup", status="done", removed=sample_scratch)
            except Exception as exc:
                failures[run_id] = str(exc)
                logger.emit(sample=run_id, step="profile-sra", status="failed", error=exc)
                logger.emit(
                    sample=run_id,
                    step="checkpoint",
                    status="retained",
                    path=Path(scratch_dir) / run_id,
                )
    logger.emit(
        step="profile-sra",
        status="done",
        samples=len(run_ids),
        completed=len(run_ids) - len(failures),
        failed=len(failures),
    )
    if failures and raise_on_error:
        details = "; ".join(f"{run_id}: {error}" for run_id, error in sorted(failures.items()))
        raise RuntimeError(f"{len(failures)} sample profiling job(s) failed: {details}")
    return failures


def sync_remaining_profiles(
    *,
    db_file: Path,
    cache_dir: Path,
    scratch_dir: Path,
    output_dir: Path,
    sylph_db: Path | None = None,
    accessions_dir: Path | None = None,
    threads: int = 8,
    workflow_config: WorkflowConfig | None = None,
    check_dependencies: bool = True,
    cleanup_outputs: bool = True,
    logger: WorkflowLogger | None = None,
) -> SyncSummary:
    """Profile remaining SRA runs and import completed output bundles."""
    logger = logger or WorkflowLogger()
    if check_dependencies:
        logger.emit(step="dependency-check", status="start")
        healthcheck.assert_dependencies()
        logger.emit(step="dependency-check", status="done")

    with db.connect(db_file) as conn:
        run_ids = db.remaining_runs(conn)

    if not run_ids:
        logger.emit(step="sync", status="done", remaining=0, imported=0)
        return SyncSummary(requested=0, imported=0, skipped=0, failed=0, cleaned_files=0)

    logger.emit(step="sync", status="start", remaining=len(run_ids), output_dir=output_dir)
    imported = 0
    skipped = 0
    cleaned_files = 0
    with db.connect(db_file) as conn:
        def import_completed_sample(run_id: str) -> None:
            nonlocal imported, cleaned_files
            try:
                bundle = discover_profile_bundle(output_dir=output_dir, run_id=run_id)
            except FileNotFoundError as exc:
                logger.emit(sample=run_id, step="import", status="missing-output", error=exc)
                raise
            logger.emit(sample=run_id, step="import", status="start")
            db.import_profile_bundle(conn, bundle)
            imported += 1
            logger.emit(sample=run_id, step="import", status="done", imported=imported)
            if cleanup_outputs:
                removed = cleanup_profile_bundle(bundle)
                cleaned_files += removed
                logger.emit(sample=run_id, step="cleanup", status="done", removed_files=removed)

        failures = profile_sra_runs(
            run_ids=run_ids,
            db_file=db_file,
            cache_dir=cache_dir,
            scratch_dir=scratch_dir,
            sylph_db=sylph_db,
            output_dir=output_dir,
            accessions_dir=accessions_dir,
            threads=threads,
            workflow_config=workflow_config,
            logger=logger,
            raise_on_error=False,
            completion_callback=import_completed_sample,
        ) or {}
        for run_id, error in failures.items():
            db.mark_run_failed(conn, run_id=run_id, error=error)

    failed = len(failures)
    logger.emit(
        step="sync",
        status="done",
        remaining=len(run_ids),
        imported=imported,
        skipped=skipped,
        failed=failed,
        cleaned_files=cleaned_files,
    )
    return SyncSummary(
        requested=len(run_ids),
        imported=imported,
        skipped=skipped,
        failed=failed,
        cleaned_files=cleaned_files,
    )


def discover_profile_bundle(*, output_dir: Path, run_id: str) -> db.ProfileBundle:
    """Find the conventional output files for one profiled SRA run."""
    output_dir = Path(output_dir)
    profile_file = _first_existing_path(
        output_dir,
        [
            f"{run_id}.profile.parquet",
            f"{run_id}.parquet",
        ],
        required_name="profile_file",
        run_id=run_id,
    )
    genome_stats_file = _first_existing_path(
        output_dir,
        [
            f"{run_id}.genome_stats.parquet",
            f"{run_id}_genome_stats.parquet",
        ],
        required_name="genome_stats_file",
        run_id=run_id,
    )
    sylph_abundance_file = _first_existing_path(
        output_dir,
        [
            f"{run_id}.sylph.parquet",
            f"{run_id}.sylph.csv",
            f"{run_id}.sylph.tsv",
            f"{run_id}_sylph.parquet",
            f"{run_id}_sylph.csv",
            f"{run_id}_sylph.tsv",
        ],
        required_name="sylph_abundance_file",
        run_id=run_id,
    )
    gene_stats_file = _optional_existing_path(
        output_dir,
        [
            f"{run_id}.gene_stats.parquet",
            f"{run_id}_gene_stats.parquet",
        ],
    )
    return db.ProfileBundle(
        run_id=run_id,
        profile_file=profile_file,
        genome_stats_file=genome_stats_file,
        gene_stats_file=gene_stats_file,
        sylph_abundance_file=sylph_abundance_file,
    )


def cleanup_profile_bundle(bundle: db.ProfileBundle) -> int:
    """Delete per-sample output files after they have been imported."""
    removed = 0
    paths = [
        bundle.profile_file,
        bundle.genome_stats_file,
        bundle.gene_stats_file,
        bundle.sylph_abundance_file,
    ]
    seen: set[Path] = set()
    for path in paths:
        if path is None:
            continue
        clean_path = Path(path)
        if clean_path in seen:
            continue
        seen.add(clean_path)
        if clean_path.exists():
            clean_path.unlink()
            removed += 1
    return removed


def _profile_one_sra_run(
    *,
    run_id: str,
    cache_manager: cache.GenomeCache,
    scratch_dir: Path,
    sylph_db: Path | None,
    output_dir: Path | None,
    accessions_dir: Path | None,
    runtime: WorkflowRuntime,
    logger: WorkflowLogger,
) -> None:
    sample_scratch = Path(scratch_dir) / run_id
    sample_scratch.mkdir(parents=True, exist_ok=True)
    reads_dir = sample_scratch / "reads"
    reads_dir.mkdir(parents=True, exist_ok=True)
    sra_dir = sample_scratch / "sra"
    sra_archive = sra_dir / run_id / f"{run_id}.sra"

    if not _sample_fastqs(sample_scratch, required=False):
        logger.emit(sample=run_id, step="download", status="start")
        download_threads = runtime.threads("sra_download")
        if not _valid_file(sra_archive):
            sra_dir.mkdir(parents=True, exist_ok=True)
            runtime.run(
                "sra_download",
                ["prefetch", "--output-directory", str(sra_dir), run_id],
                sample=run_id,
            )
        else:
            logger.emit(sample=run_id, step="prefetch", status="cached", file=sra_archive)
        runtime.run(
            "sra_download",
            ["fasterq-dump", str(sra_archive), "--outdir", str(reads_dir), "--threads", str(download_threads)],
            sample=run_id,
        )
        logger.emit(sample=run_id, step="download", status="done", fastqs=len(_sample_fastqs(sample_scratch)))
    else:
        logger.emit(sample=run_id, step="download", status="cached", fastqs=len(_sample_fastqs(sample_scratch)))

    accessions_file = sample_scratch / "accessions.txt"
    if not accessions_file.exists() and accessions_dir is not None:
        source = _find_sample_accessions_file(accessions_dir, run_id)
        if source is not None:
            shutil.copy2(source, accessions_file)
            logger.emit(sample=run_id, step="sylph", status="using-accessions-file", file=source)
    sylph_output = sample_scratch / f"{run_id}.sylph.tsv"
    if not accessions_file.exists():
        if sylph_db is None:
            raise FileNotFoundError(
                f"sample={run_id} step=sylph missing accession list and no Sylph database was provided. "
                "Pass --sylph-db so MetaTrawl can run `sylph profile`, or pass --accessions-dir as a manual override."
            )
        logger.emit(sample=run_id, step="sylph", status="start", syldb=sylph_db)
        _run_sylph_profile(
            sylph_db=sylph_db,
            sample_scratch=sample_scratch,
            run_id=run_id,
            runtime=runtime,
            output_file=sylph_output,
        )
        accessions = extract_accessions_from_sylph_table(sylph_output)
        if not accessions:
            raise ValueError(f"sample={run_id} step=sylph produced no nonzero genome accessions: {sylph_output}")
        _write_accessions_file(accessions_file, accessions)
        logger.emit(sample=run_id, step="sylph", status="done", genomes=len(accessions), output=sylph_output)
    else:
        logger.emit(sample=run_id, step="sylph", status="cached", file=accessions_file)
    if output_dir is not None and _valid_file(sylph_output):
        output_dir.mkdir(parents=True, exist_ok=True)
        published_sylph = output_dir / f"{run_id}.sylph.tsv"
        if not _valid_file(published_sylph):
            shutil.copy2(sylph_output, published_sylph)
            logger.emit(sample=run_id, step="sylph", status="published", file=published_sylph)

    accessions = cache.read_accessions_file(accessions_file)
    reference_dir = sample_scratch / "reference"
    reference = _cached_reference(reference_dir)
    if reference is None:
        reference = cache_manager.prepare_reference(
            accessions=accessions,
            output_dir=reference_dir,
            sample=run_id,
        )
    else:
        logger.emit(sample=run_id, step="cache", status="cached-reference", reference=reference.reference_fasta)
    logger.emit(sample=run_id, step="profile", status="ready", reference=reference.reference_fasta)
    if output_dir is not None:
        _run_alignment_and_profile(
            run_id=run_id,
            sample_scratch=sample_scratch,
            reference=reference,
            output_dir=output_dir,
            runtime=runtime,
            logger=logger,
        )


def _valid_file(path: Path) -> bool:
    return path.is_file() and path.stat().st_size > 0


def _cached_reference(reference_dir: Path) -> cache.PreparedReference | None:
    reference = cache.PreparedReference(
        reference_fasta=reference_dir / "reference.fna",
        gene_fasta=reference_dir / "genes.fna",
        stb_file=reference_dir / "reference.stb",
    )
    if all(_valid_file(path) for path in (reference.reference_fasta, reference.gene_fasta, reference.stb_file)):
        return reference
    return None


def _run(cmd: list[str], *, sample: str, step: str) -> None:
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as exc:
        raise RuntimeError(
            f"sample={sample} step={step} command failed: {' '.join(cmd)}\n{exc.stderr}"
        ) from exc


class _DirectRuntime:
    """Compatibility adapter for direct/private helper use and focused tests."""

    def __init__(self, threads: int) -> None:
        self._threads = threads

    def threads(self, stage: str) -> int:
        return self._threads

    def run(self, stage: str, cmd: list[str], *, sample: str, stdout_file: Path | None = None) -> None:
        if stdout_file is None:
            _run(cmd, sample=sample, step=stage.replace("_", "-"))
            return
        try:
            with stdout_file.open("w") as handle:
                subprocess.run(cmd, check=True, stdout=handle, stderr=subprocess.PIPE, text=True)
        except subprocess.CalledProcessError as exc:
            raise RuntimeError(f"sample={sample} step={stage} command failed: {' '.join(cmd)}\n{exc.stderr}") from exc

    def run_shell(self, stage: str, command: str, *, sample: str) -> None:
        _run_alignment_shell(command=command, sample=sample)


def _run_alignment_and_profile(
    *,
    run_id: str,
    sample_scratch: Path,
    reference: cache.PreparedReference,
    output_dir: Path,
    runtime: WorkflowRuntime | _DirectRuntime | None = None,
    threads: int | None = None,
    logger: WorkflowLogger,
) -> None:
    runtime = runtime or _DirectRuntime(threads or 1)
    if reference.stb_file is None:
        raise ValueError(f"sample={run_id} step=profile missing STB file from cache reference preparation")
    output_dir.mkdir(parents=True, exist_ok=True)
    profile_work = sample_scratch / "zipstrain_profile"
    profile_work.mkdir(parents=True, exist_ok=True)

    prepared_outputs = (
        profile_work / "genomes_bed_file.bed",
        profile_work / "gene_range_table.tsv",
        profile_work / "profiling_contract.json",
    )
    if not all(_valid_file(path) for path in prepared_outputs):
        logger.emit(sample=run_id, step="prepare-profile", status="start")
        runtime.run(
            "prepare_profile",
            [
                "zipstrain",
                "utilities",
                "prepare_profiling",
                "--reference-fasta",
                str(reference.reference_fasta),
                "--gene-fasta",
                str(reference.gene_fasta),
                "--stb-file",
                str(reference.stb_file),
                "--output-dir",
                str(profile_work),
            ],
            sample=run_id,
        )
        logger.emit(sample=run_id, step="prepare-profile", status="done")
    else:
        logger.emit(sample=run_id, step="prepare-profile", status="cached")
    _ensure_null_model(profile_work=profile_work, run_id=run_id, logger=logger, runtime=runtime)
    _ensure_reference_fai(profile_work=profile_work, run_id=run_id, logger=logger, runtime=runtime)

    if not _bowtie_index_complete(reference.reference_fasta):
        logger.emit(sample=run_id, step="bowtie2-build", status="start")
        bowtie_threads = runtime.threads("bowtie_build")
        runtime.run(
            "bowtie_build",
            ["bowtie2-build", "--threads", str(bowtie_threads), str(reference.reference_fasta), str(reference.reference_fasta)],
            sample=run_id,
        )
        logger.emit(sample=run_id, step="bowtie2-build", status="done")
    else:
        logger.emit(sample=run_id, step="bowtie2-build", status="cached")

    bam_file = sample_scratch / f"{run_id}.bam"
    if not _valid_file(bam_file):
        logger.emit(sample=run_id, step="align", status="start", bam=bam_file)
        temporary_bam = bam_file.with_suffix(".bam.tmp")
        temporary_bam.unlink(missing_ok=True)
        runtime.run_shell(
            "alignment",
            _build_alignment_command(
                reference_fasta=reference.reference_fasta,
                fastqs=_sample_fastqs(sample_scratch),
                threads=runtime.threads("alignment"),
                bam_file=temporary_bam,
            ),
            sample=run_id,
        )
        if not _valid_file(temporary_bam):
            raise RuntimeError(f"sample={run_id} step=align did not produce a non-empty BAM: {temporary_bam}")
        temporary_bam.replace(bam_file)
        logger.emit(sample=run_id, step="align", status="done", bam=bam_file)
    else:
        logger.emit(sample=run_id, step="align", status="cached", bam=bam_file)

    if not _published_profile_bundle_complete(run_id=run_id, output_dir=output_dir):
        logger.emit(sample=run_id, step="profile", status="start")
        runtime.run(
            "profile",
            [
                "zipstrain",
                "utilities",
                "profile-single",
                "--bam-file",
                str(bam_file),
                "--bed-file",
                str(profile_work / "genomes_bed_file.bed"),
                "--stb-file",
                str(reference.stb_file),
                "--null-model",
                str(profile_work / "null_model.parquet"),
                "--gene-range-table",
                str(profile_work / "gene_range_table.tsv"),
                "--reference-fasta",
                str(profile_work / "reference.fasta"),
                "--profiling-contract",
                str(profile_work / "profiling_contract.json"),
                "--max-concurrency",
                str(runtime.threads("profile")),
                "--output-dir",
                str(profile_work),
            ],
            sample=run_id,
        )
        _publish_profile_outputs(run_id=run_id, profile_work=profile_work, output_dir=output_dir)
        logger.emit(sample=run_id, step="profile", status="done", output_dir=output_dir)
    else:
        logger.emit(sample=run_id, step="profile", status="cached", output_dir=output_dir)


def _sample_fastqs(sample_scratch: Path, *, required: bool = True) -> list[Path]:
    fastqs: list[Path] = []
    for pattern in ("*.fastq", "*.fq", "*.fastq.gz", "*.fq.gz"):
        fastqs.extend(sorted(sample_scratch.rglob(pattern)))
    fastqs = [path for path in fastqs if _valid_file(path)]
    if required and not fastqs:
        raise FileNotFoundError(f"step=align found no FASTQ files in {sample_scratch}")
    return fastqs


def _bowtie_index_complete(reference_fasta: Path) -> bool:
    small = [Path(f"{reference_fasta}.{suffix}.bt2") for suffix in ("1", "2", "3", "4", "rev.1", "rev.2")]
    large = [Path(f"{reference_fasta}.{suffix}.bt2l") for suffix in ("1", "2", "3", "4", "rev.1", "rev.2")]
    return all(_valid_file(path) for path in small) or all(_valid_file(path) for path in large)


def _published_profile_bundle_complete(*, run_id: str, output_dir: Path) -> bool:
    required = (
        output_dir / f"{run_id}.profile.parquet",
        output_dir / f"{run_id}.genome_stats.parquet",
        output_dir / f"{run_id}.gene_stats.parquet",
        output_dir / f"{run_id}.sylph.tsv",
    )
    return all(_valid_file(path) for path in required)


def _ensure_null_model(*, profile_work: Path, run_id: str, logger: WorkflowLogger, runtime: WorkflowRuntime) -> Path:
    null_model = profile_work / "null_model.parquet"
    if null_model.exists():
        return null_model
    logger.emit(sample=run_id, step="null-model", status="start", output=null_model)
    runtime.run(
        "prepare_profile",
        [
            "zipstrain",
            "utilities",
            "build-null-model",
            "--output-file",
            str(null_model),
        ],
        sample=run_id,
    )
    logger.emit(sample=run_id, step="null-model", status="done", output=null_model)
    return null_model


def _ensure_reference_fai(*, profile_work: Path, run_id: str, logger: WorkflowLogger, runtime: WorkflowRuntime) -> Path:
    """Build the FASTA index once before concurrent mpileup workers can race on it."""
    reference_fasta = profile_work / "reference.fasta"
    if not reference_fasta.exists():
        raise FileNotFoundError(f"sample={run_id} step=reference-index missing reference FASTA: {reference_fasta}")
    fai_file = Path(str(reference_fasta) + ".fai")
    logger.emit(sample=run_id, step="reference-index", status="start", reference=reference_fasta)
    fai_file.unlink(missing_ok=True)
    runtime.run(
        "prepare_profile",
        ["samtools", "faidx", str(reference_fasta)],
        sample=run_id,
    )
    if not _valid_file(fai_file):
        raise RuntimeError(f"sample={run_id} step=reference-index did not produce FASTA index: {fai_file}")
    logger.emit(sample=run_id, step="reference-index", status="done", index=fai_file)
    return fai_file


def _build_alignment_command(*, reference_fasta: Path, fastqs: list[Path], threads: int, bam_file: Path) -> str:
    quoted_reference = shlex.quote(str(reference_fasta))
    quoted_threads = shlex.quote(str(threads))
    if len(fastqs) == 2:
        read_args = f"-1 {shlex.quote(str(fastqs[0]))} -2 {shlex.quote(str(fastqs[1]))}"
    else:
        read_args = f"-U {shlex.quote(','.join(str(path) for path in fastqs))}"
    return (
        f"bowtie2 -x {quoted_reference} {read_args} --threads {quoted_threads} "
        f"| samtools view -bS -F 4 - "
        f"| samtools sort -@ {quoted_threads} -o {shlex.quote(str(bam_file))} -"
    )


def _run_alignment_shell(*, command: str, sample: str) -> None:
    try:
        subprocess.run(command, check=True, shell=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as exc:
        raise RuntimeError(f"sample={sample} step=align command failed: {command}\n{exc.stderr}") from exc


def _publish_profile_outputs(*, run_id: str, profile_work: Path, output_dir: Path) -> None:
    outputs = {
        profile_work / f"{run_id}_profile.parquet": output_dir / f"{run_id}.profile.parquet",
        profile_work / f"{run_id}_genome_stats.parquet": output_dir / f"{run_id}.genome_stats.parquet",
        profile_work / f"{run_id}_gene_stats.parquet": output_dir / f"{run_id}.gene_stats.parquet",
    }
    missing = [str(source) for source in outputs if not source.exists()]
    if missing:
        raise FileNotFoundError(f"sample={run_id} step=profile missing expected ZipStrain outputs: {', '.join(missing)}")
    for source, destination in outputs.items():
        shutil.copy2(source, destination)


def _run_sylph_profile(
    *, sylph_db: Path, sample_scratch: Path, run_id: str,
    runtime: WorkflowRuntime | _DirectRuntime | None = None, threads: int | None = None, output_file: Path
) -> None:
    runtime = runtime or _DirectRuntime(threads or 1)
    fastqs = _sample_fastqs(sample_scratch)
    if len(fastqs) == 2:
        cmd = ["sylph", "profile", str(sylph_db), "-1", str(fastqs[0]), "-2", str(fastqs[1]), "-t", str(runtime.threads("sylph"))]
    else:
        cmd = ["sylph", "profile", str(sylph_db), *[str(path) for path in fastqs], "-t", str(runtime.threads("sylph"))]
    runtime.run("sylph", cmd, sample=run_id, stdout_file=output_file)


def _resolve_existing_sylph_db(sylph_db: Path | str) -> Path:
    path = Path(sylph_db).expanduser()
    if not path.exists():
        raise FileNotFoundError(
            f"step=sylph Sylph database does not exist: {path}. "
            "Pass an absolute --sylph-db path or run from the directory containing the .syldb file."
        )
    return path.resolve()


def extract_accessions_from_sylph_table(sylph_output: Path) -> list[str]:
    table = _read_sylph_table(sylph_output)
    if table.is_empty():
        return []
    genome_col = _first_matching_column(table.columns, SYLPH_GENOME_COLUMNS)
    if genome_col is None:
        raise ValueError(
            f"step=sylph could not find a genome column in {sylph_output}. "
            f"Expected one of: {', '.join(SYLPH_GENOME_COLUMNS)}"
        )
    table = _filter_nonzero_abundance_rows(table)
    if table.is_empty():
        return []
    accessions = (
        table.select(
            pl.col(genome_col)
            .cast(pl.Utf8)
            .str.extract(ACCESSION_PATTERN.pattern, 1)
            .str.to_uppercase()
            .alias("accession")
        )
        .drop_nulls("accession")
        .unique()
        .sort("accession")
        .get_column("accession")
        .to_list()
    )
    return [str(accession) for accession in accessions]


def _read_sylph_table(path: Path) -> pl.DataFrame:
    suffix = path.suffix.lower()
    if suffix == ".parquet":
        return pl.read_parquet(path)
    if suffix in {".tsv", ".txt"}:
        return pl.read_csv(path, separator="\t")
    return pl.read_csv(path)


def _filter_nonzero_abundance_rows(table: pl.DataFrame) -> pl.DataFrame:
    numeric_cols = [
        name
        for name, dtype in table.schema.items()
        if dtype.is_numeric()
    ]
    if not numeric_cols:
        return table
    preferred = _first_matching_column(numeric_cols, SYLPH_ABUNDANCE_COLUMNS)
    if preferred is not None:
        cols_to_check = [preferred]
    else:
        abundance_like = [name for name in numeric_cols if ("abund" in name.lower()) or ("cov" in name.lower())]
        cols_to_check = abundance_like if abundance_like else numeric_cols
    return table.filter(
        pl.any_horizontal([pl.col(name).cast(pl.Float64, strict=False).fill_null(0.0) > 0.0 for name in cols_to_check])
    )


def _write_accessions_file(path: Path, accessions: list[str]) -> None:
    path.write_text("".join(f"{accession}\n" for accession in accessions))


def _first_matching_column(columns: list[str], preferred: list[str]) -> str | None:
    by_lower = {column.lower(): column for column in columns}
    for name in preferred:
        match = by_lower.get(name.lower())
        if match is not None:
            return match
    return None


def _find_sample_accessions_file(accessions_dir: Path, run_id: str) -> Path | None:
    return _optional_existing_path(
        Path(accessions_dir),
        [
            f"{run_id}.accessions.txt",
            f"{run_id}.accessions.csv",
            f"{run_id}.txt",
            f"{run_id}.csv",
            f"{run_id}/accessions.txt",
            f"{run_id}/accessions.csv",
        ],
    )


def _first_existing_path(output_dir: Path, names: list[str], *, required_name: str, run_id: str) -> Path:
    found = _optional_existing_path(output_dir, names)
    if found is None:
        candidates = ", ".join(str(output_dir / name) for name in names)
        raise FileNotFoundError(f"sample={run_id} step=import missing {required_name}; looked for: {candidates}")
    return found


def _optional_existing_path(output_dir: Path, names: list[str]) -> Path | None:
    for name in names:
        path = output_dir / name
        if path.exists():
            return path
    return None

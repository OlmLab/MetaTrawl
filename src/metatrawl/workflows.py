"""Workflow adapters for MetaTrawl."""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
import re
import shlex
import shutil
import subprocess
from tempfile import TemporaryDirectory
from dataclasses import dataclass

import polars as pl

from metatrawl import cache
from metatrawl import db
from metatrawl import healthcheck
from metatrawl.logging import WorkflowLogger


ACCESSION_PATTERN = re.compile(r"\b((?:GC[AF])_\d+(?:\.\d+)?)", re.IGNORECASE)
SYLPH_GENOME_COLUMNS = ["Genome_file", "genome_file", "genome", "Genome", "Reference", "reference", "target"]
SYLPH_ABUNDANCE_COLUMNS = ["eff_cov", "genome_cov", "abundance", "relative_abundance", "rel_abundance", "Taxonomic_abundance"]


@dataclass(frozen=True)
class SyncSummary:
    """High-level result for one MetaTrawl sync run."""

    requested: int
    imported: int
    skipped: int
    cleaned_files: int


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
    if db.get_matrix_store(conn, chosen_matrix_id) is not None and not overwrite:
        raise ValueError(f"Matrix ID already exists: {chosen_matrix_id}")

    from zipstrain import matrix_pairs as mp

    logger.emit(step="matrix-build", status="exporting-profiles", samples=len(sample_ids), genome=genome)
    with TemporaryDirectory(prefix="metatrawl_matrix_profiles_") as tmp_dir:
        profile_dir = Path(tmp_dir)
        db.export_profile_parquets(conn, sample_ids=sample_ids, output_dir=profile_dir)
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
        overwrite=overwrite,
    )
    logger.emit(step="matrix-build", status="done", matrix_id=chosen_matrix_id, samples=len(sample_ids))
    return store


def append_matrix_from_database(
    conn,
    *,
    matrix_id: str,
    memory_limit_gb: float = 16.0,
    export_batch_mb: float = 128.0,
    logger: WorkflowLogger | None = None,
) -> int:
    """Append newly imported complete samples into a registered matrix store."""
    logger = logger or WorkflowLogger()
    store = db.get_matrix_store(conn, matrix_id)
    if store is None:
        raise ValueError(f"Unknown matrix ID: {matrix_id}")
    sample_ids = db.unmaterialized_sample_ids(conn, matrix_id)
    if not sample_ids:
        raise ValueError(f"No new complete samples are available to append to matrix ID: {matrix_id}")

    from zipstrain import matrix_pairs as mp

    logger.emit(step="matrix-append", status="exporting-profiles", matrix_id=matrix_id, samples=len(sample_ids))
    with TemporaryDirectory(prefix="metatrawl_matrix_append_") as tmp_dir:
        profile_dir = Path(tmp_dir)
        db.export_profile_parquets(conn, sample_ids=sample_ids, output_dir=profile_dir)
        logger.emit(step="matrix-append", status="appending", matrix_id=matrix_id, samples=len(sample_ids))
        mp.append_matrix_hdf5(
            profile_dir=profile_dir,
            matrix_hdf5_file=store.matrix_file,
            memory_limit_gb=memory_limit_gb,
            export_batch_mb=export_batch_mb,
        )
    db.add_matrix_store_samples(conn, matrix_id=matrix_id, sample_ids=sample_ids)
    db.update_matrix_profile_count(conn, matrix_id=matrix_id)
    logger.emit(step="matrix-append", status="done", matrix_id=matrix_id, samples=len(sample_ids))
    return len(sample_ids)


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
    logger: WorkflowLogger | None = None,
) -> None:
    """Download/profile SRA runs, import outputs, and delete per-sample scratch.

    This is intentionally conservative: it wires the lifecycle and cleanup, while
    external-command behavior remains visible and easy to replace in tests.
    """
    logger = logger or WorkflowLogger()
    cache_manager = cache.GenomeCache(cache_dir, logger=logger)
    sylph_db = _resolve_existing_sylph_db(sylph_db) if sylph_db is not None else None
    max_workers = max(1, min(len(run_ids), threads))
    logger.emit(step="profile-sra", status="start", samples=len(run_ids), workers=max_workers, db=db_file)
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = [
            executor.submit(
                _profile_one_sra_run,
                run_id=run_id,
                cache_manager=cache_manager,
                scratch_dir=Path(scratch_dir),
                sylph_db=sylph_db,
                output_dir=Path(output_dir) if output_dir is not None else None,
                accessions_dir=Path(accessions_dir) if accessions_dir is not None else None,
                threads=max(1, threads // max_workers),
                logger=logger,
            )
            for run_id in run_ids
        ]
        for future in as_completed(futures):
            future.result()
    logger.emit(step="profile-sra", status="done", samples=len(run_ids))


def sync_remaining_profiles(
    *,
    db_file: Path,
    cache_dir: Path,
    scratch_dir: Path,
    output_dir: Path,
    sylph_db: Path | None = None,
    accessions_dir: Path | None = None,
    threads: int = 8,
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
        return SyncSummary(requested=0, imported=0, skipped=0, cleaned_files=0)

    logger.emit(step="sync", status="start", remaining=len(run_ids), output_dir=output_dir)
    profile_sra_runs(
        run_ids=run_ids,
        db_file=db_file,
        cache_dir=cache_dir,
        scratch_dir=scratch_dir,
        sylph_db=sylph_db,
        output_dir=output_dir,
        accessions_dir=accessions_dir,
        threads=threads,
        logger=logger,
    )

    imported = 0
    skipped = 0
    cleaned_files = 0
    with db.connect(db_file) as conn:
        for run_id in run_ids:
            try:
                bundle = discover_profile_bundle(output_dir=output_dir, run_id=run_id)
            except FileNotFoundError as exc:
                skipped += 1
                logger.emit(sample=run_id, step="import", status="missing-output", error=exc)
                raise
            logger.emit(sample=run_id, step="import", status="start")
            db.import_profile_bundle(conn, bundle)
            imported += 1
            logger.emit(sample=run_id, step="import", status="done")
            if cleanup_outputs:
                removed = cleanup_profile_bundle(bundle)
                cleaned_files += removed
                logger.emit(sample=run_id, step="cleanup", status="done", removed_files=removed)

    logger.emit(step="sync", status="done", remaining=len(run_ids), imported=imported, skipped=skipped, cleaned_files=cleaned_files)
    return SyncSummary(requested=len(run_ids), imported=imported, skipped=skipped, cleaned_files=cleaned_files)


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
    threads: int,
    logger: WorkflowLogger,
) -> None:
    sample_scratch = Path(scratch_dir) / run_id
    sample_scratch.mkdir(parents=True, exist_ok=True)
    try:
        logger.emit(sample=run_id, step="download", status="start")
        _run(["prefetch", run_id], sample=run_id, step="download")
        _run(["fasterq-dump", run_id, "--outdir", str(sample_scratch), "--threads", str(threads)], sample=run_id, step="download")
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
                threads=threads,
                output_file=sylph_output,
            )
            accessions = extract_accessions_from_sylph_table(sylph_output)
            if not accessions:
                raise ValueError(f"sample={run_id} step=sylph produced no nonzero genome accessions: {sylph_output}")
            _write_accessions_file(accessions_file, accessions)
            logger.emit(sample=run_id, step="sylph", status="done", genomes=len(accessions), output=sylph_output)
            if output_dir is not None:
                output_dir.mkdir(parents=True, exist_ok=True)
                published_sylph = output_dir / f"{run_id}.sylph.tsv"
                shutil.copy2(sylph_output, published_sylph)
                logger.emit(sample=run_id, step="sylph", status="published", file=published_sylph)
        if not accessions_file.exists():
            raise FileNotFoundError(
                f"sample={run_id} step=sylph missing accession list. "
                f"Expected scratch file: {accessions_file}"
            )
        accessions = cache.read_accessions_file(accessions_file)
        reference = cache_manager.prepare_reference(
            accessions=accessions,
            output_dir=sample_scratch / "reference",
            sample=run_id,
        )
        logger.emit(sample=run_id, step="profile", status="ready", reference=reference.reference_fasta)
        if output_dir is not None:
            _run_alignment_and_profile(
                run_id=run_id,
                sample_scratch=sample_scratch,
                reference=reference,
                output_dir=output_dir,
                threads=threads,
                logger=logger,
            )
    finally:
        if sample_scratch.exists():
            shutil.rmtree(sample_scratch)
            logger.emit(sample=run_id, step="cleanup", status="done", removed=sample_scratch)


def _run(cmd: list[str], *, sample: str, step: str) -> None:
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as exc:
        raise RuntimeError(
            f"sample={sample} step={step} command failed: {' '.join(cmd)}\n{exc.stderr}"
        ) from exc


def _run_alignment_and_profile(
    *,
    run_id: str,
    sample_scratch: Path,
    reference: cache.PreparedReference,
    output_dir: Path,
    threads: int,
    logger: WorkflowLogger,
) -> None:
    if reference.stb_file is None:
        raise ValueError(f"sample={run_id} step=profile missing STB file from cache reference preparation")
    output_dir.mkdir(parents=True, exist_ok=True)
    profile_work = sample_scratch / "zipstrain_profile"
    profile_work.mkdir(parents=True, exist_ok=True)

    logger.emit(sample=run_id, step="prepare-profile", status="start")
    _run(
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
        step="prepare-profile",
    )
    logger.emit(sample=run_id, step="prepare-profile", status="done")

    logger.emit(sample=run_id, step="bowtie2-build", status="start")
    _run(["bowtie2-build", "--threads", str(threads), str(reference.reference_fasta), str(reference.reference_fasta)], sample=run_id, step="bowtie2-build")
    logger.emit(sample=run_id, step="bowtie2-build", status="done")

    bam_file = sample_scratch / f"{run_id}.bam"
    logger.emit(sample=run_id, step="align", status="start", bam=bam_file)
    _run_alignment_shell(
        command=_build_alignment_command(
            reference_fasta=reference.reference_fasta,
            fastqs=_sample_fastqs(sample_scratch),
            threads=threads,
            bam_file=bam_file,
        ),
        sample=run_id,
    )
    logger.emit(sample=run_id, step="align", status="done", bam=bam_file)

    logger.emit(sample=run_id, step="profile", status="start")
    _run(
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
            "--max-concurrency",
            str(threads),
            "--output-dir",
            str(profile_work),
        ],
        sample=run_id,
        step="profile",
    )
    _publish_profile_outputs(run_id=run_id, profile_work=profile_work, output_dir=output_dir)
    logger.emit(sample=run_id, step="profile", status="done", output_dir=output_dir)


def _sample_fastqs(sample_scratch: Path) -> list[Path]:
    fastqs = sorted(sample_scratch.glob("*.fastq")) + sorted(sample_scratch.glob("*.fq"))
    fastqs += sorted(sample_scratch.glob("*.fastq.gz")) + sorted(sample_scratch.glob("*.fq.gz"))
    if not fastqs:
        raise FileNotFoundError(f"step=align found no FASTQ files in {sample_scratch}")
    return fastqs


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


def _run_sylph_profile(*, sylph_db: Path, sample_scratch: Path, run_id: str, threads: int, output_file: Path) -> None:
    fastqs = sorted(sample_scratch.glob("*.fastq")) + sorted(sample_scratch.glob("*.fq"))
    fastqs += sorted(sample_scratch.glob("*.fastq.gz")) + sorted(sample_scratch.glob("*.fq.gz"))
    if not fastqs:
        raise FileNotFoundError(f"sample={run_id} step=sylph found no FASTQ files in {sample_scratch}")
    if len(fastqs) == 2:
        cmd = ["sylph", "profile", str(sylph_db), "-1", str(fastqs[0]), "-2", str(fastqs[1]), "-t", str(threads)]
    else:
        cmd = ["sylph", "profile", str(sylph_db), *[str(path) for path in fastqs], "-t", str(threads)]
    try:
        with output_file.open("w") as handle:
            subprocess.run(cmd, check=True, stdout=handle, stderr=subprocess.PIPE, text=True)
    except subprocess.CalledProcessError as exc:
        raise RuntimeError(f"sample={run_id} step=sylph command failed: {' '.join(cmd)}\n{exc.stderr}") from exc


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

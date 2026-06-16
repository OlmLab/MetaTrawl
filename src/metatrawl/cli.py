"""MetaTrawl command-line interface."""

from __future__ import annotations

import csv
import json
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

import click
from rich.console import Console
from rich.table import Table

from metatrawl import __version__
from metatrawl import cache
from metatrawl import db as registry
from metatrawl import healthcheck
from metatrawl import workflows
from metatrawl.logging import WorkflowLogger


@click.group()
@click.version_option(version=__version__, prog_name="metatrawl")
def cli() -> None:
    """Registry tooling for SRA-scale ZipStrain profiling projects."""


@cli.command("init")
@click.option("--db", "db_file", required=True, type=click.Path(path_type=Path), help="MetaTrawl DuckDB registry.")
def init(db_file: Path) -> None:
    """Initialize a MetaTrawl project database."""
    conn = registry.connect(db_file)
    conn.close()
    click.echo(f"initialized={db_file}")


@cli.group("runs")
def runs_group() -> None:
    """Manage registered SRA run IDs."""


@runs_group.command("add")
@click.option("--db", "db_file", required=True, type=click.Path(path_type=Path), help="MetaTrawl DuckDB registry.")
@click.argument("run_ids", nargs=-1, required=True)
def runs_add(db_file: Path, run_ids: tuple[str, ...]) -> None:
    """Add SRA run IDs to the registry."""
    with registry.connect(db_file) as conn:
        added, reactivated = registry.add_runs(conn, list(run_ids))
    click.echo(f"added={added} reactivated={reactivated}")


@runs_group.command("list")
@click.option("--db", "db_file", required=True, type=click.Path(path_type=Path), help="MetaTrawl DuckDB registry.")
@click.option("--include-deleted", is_flag=True, help="Show soft-deleted runs too.")
def runs_list(db_file: Path, include_deleted: bool) -> None:
    """List registered SRA runs."""
    with registry.connect(db_file) as conn:
        rows = registry.list_runs(conn, include_deleted=include_deleted)
    _print_rows(rows, columns=["run_id", "status", "deleted_at"], title="SRA runs")


@runs_group.command("delete")
@click.option("--db", "db_file", required=True, type=click.Path(path_type=Path), help="MetaTrawl DuckDB registry.")
@click.argument("run_ids", nargs=-1, required=True)
def runs_delete(db_file: Path, run_ids: tuple[str, ...]) -> None:
    """Soft-delete SRA run IDs from future work."""
    with registry.connect(db_file) as conn:
        deleted = registry.delete_runs(conn, list(run_ids))
    click.echo(f"deleted={deleted}")


@cli.group("profiles")
def profiles_group() -> None:
    """Manage completed ZipStrain/Sylph profile imports."""


@profiles_group.command("remaining")
@click.option("--db", "db_file", required=True, type=click.Path(path_type=Path), help="MetaTrawl DuckDB registry.")
@click.option("--output-file", "-o", type=click.Path(path_type=Path), help="CSV path for remaining SRA run IDs.")
def profiles_remaining(db_file: Path, output_file: Path | None) -> None:
    """Write active SRA run IDs that still need profile imports."""
    with registry.connect(db_file) as conn:
        run_ids = registry.remaining_runs(conn)
    if output_file is not None:
        _write_remaining_csv(output_file, run_ids)
        click.echo(f"wrote={output_file} remaining={len(run_ids)}")
    if not run_ids:
        click.echo("All added runs have complete profiles")
    elif output_file is None:
        for run_id in run_ids:
            click.echo(run_id)


@profiles_group.command("import")
@click.option("--db", "db_file", required=True, type=click.Path(path_type=Path), help="MetaTrawl DuckDB registry.")
@click.option("--run-id", required=True, help="SRA run ID / sample ID.")
@click.option("--profile-file", required=True, type=click.Path(path_type=Path), help="ZipStrain profile parquet.")
@click.option("--genome-stats-file", required=True, type=click.Path(path_type=Path), help="ZipStrain genome stats table.")
@click.option("--gene-stats-file", type=click.Path(path_type=Path), help="Optional ZipStrain gene stats table.")
@click.option("--sylph-abundance-file", required=True, type=click.Path(path_type=Path), help="Sylph abundance table.")
@click.option("--add-run", is_flag=True, help="Register run ID before importing if it is missing.")
def profiles_import(
    db_file: Path,
    run_id: str,
    profile_file: Path,
    genome_stats_file: Path,
    gene_stats_file: Path | None,
    sylph_abundance_file: Path,
    add_run: bool,
) -> None:
    """Import profile rows, stats, and Sylph abundance into DuckDB."""
    bundle = registry.ProfileBundle(
        run_id=run_id,
        profile_file=profile_file,
        genome_stats_file=genome_stats_file,
        gene_stats_file=gene_stats_file,
        sylph_abundance_file=sylph_abundance_file,
    )
    try:
        with registry.connect(db_file) as conn:
            registry.import_profile_bundle(conn, bundle, add_run_if_missing=add_run)
    except (FileNotFoundError, ValueError) as exc:
        raise click.ClickException(str(exc)) from exc
    click.echo(f"imported={run_id}")


@profiles_group.command("add")
@click.option("--db", "db_file", required=True, type=click.Path(path_type=Path), help="MetaTrawl DuckDB registry.")
@click.option("--manifest", required=True, type=click.Path(path_type=Path), help="Completed profile manifest CSV.")
@click.option("--add-runs", is_flag=True, help="Register unknown run IDs before adding profile bundles.")
def profiles_add(db_file: Path, manifest: Path, add_runs: bool) -> None:
    """Import completed profile bundles from a manifest CSV."""
    bundles = _read_profile_manifest(manifest)
    try:
        with registry.connect(db_file) as conn:
            imported = registry.add_profiles(conn, bundles, add_runs_if_missing=add_runs)
    except (FileNotFoundError, ValueError) as exc:
        raise click.ClickException(str(exc)) from exc
    click.echo(f"imported={imported}")


@profiles_group.command("list")
@click.option("--db", "db_file", required=True, type=click.Path(path_type=Path), help="MetaTrawl DuckDB registry.")
def profiles_list(db_file: Path) -> None:
    """List imported profile bundles."""
    with registry.connect(db_file) as conn:
        rows = registry.list_profiles(conn)
    _print_rows(
        rows,
        columns=["run_id", "profile_file", "genome_stats_file", "gene_stats_file", "sylph_abundance_file"],
        title="Profiles",
    )


@cli.group("cache")
def cache_group() -> None:
    """Prepare shared genome/prodigal cache entries."""


@cache_group.command("prepare")
@click.option("--cache-dir", required=True, type=click.Path(path_type=Path), help="Shared genome cache directory.")
@click.option("--accessions", "accessions_file", required=True, type=click.Path(path_type=Path), help="Accessions text/CSV file.")
@click.option("--output-dir", required=True, type=click.Path(path_type=Path), help="Directory for temporary concatenated FASTAs.")
def cache_prepare(cache_dir: Path, accessions_file: Path, output_dir: Path) -> None:
    """Prepare temporary concatenated genome/gene FASTAs for one sample."""
    logger = WorkflowLogger()
    accessions = cache.read_accessions_file(accessions_file)
    prepared = cache.prepare_cache_reference(cache_dir=cache_dir, accessions=accessions, output_dir=output_dir, logger=logger)
    click.echo(json.dumps({"reference_fasta": str(prepared.reference_fasta), "gene_fasta": str(prepared.gene_fasta)}))


@cache_group.command("build-matrix-files")
@click.option("--genome-dir", required=True, type=click.Path(path_type=Path), help="Directory containing ACCESSION.fna files.")
@click.option("--gene-dir", required=True, type=click.Path(path_type=Path), help="Directory containing ACCESSION.genes.fna files.")
@click.option("--output-dir", required=True, type=click.Path(path_type=Path), help="Directory for generated matrix reference files.")
@click.option("--genome", help="Build files for one genome accession.")
@click.option("--accessions", "accessions_file", type=click.Path(path_type=Path), help="Optional accession list; defaults to every cached genome.")
def cache_build_matrix_files(
    genome_dir: Path,
    gene_dir: Path,
    output_dir: Path,
    genome: str | None,
    accessions_file: Path | None,
) -> None:
    """Build reusable BED, STB, and gene-range files from cache directories."""
    if genome is not None and accessions_file is not None:
        raise click.UsageError("Use either --genome or --accessions, not both.")
    accessions = cache.read_accessions_file(accessions_file) if accessions_file is not None else None
    try:
        files = cache.build_matrix_reference_files(
            genome_dir=genome_dir,
            gene_dir=gene_dir,
            output_dir=output_dir,
            accessions=accessions,
            genome=genome,
        )
    except (FileNotFoundError, ValueError) as exc:
        raise click.ClickException(str(exc)) from exc
    click.echo(
        json.dumps(
            {
                "accessions": len(files.accessions),
                "reference_fasta": str(files.reference_fasta),
                "gene_fasta": str(files.gene_fasta),
                "bed_file": str(files.bed_file),
                "stb_file": str(files.stb_file),
                "gene_range_table": str(files.gene_range_table),
            }
        )
    )


@cache_group.command("sync-matrix-files")
@click.option("--cache-dir", required=True, type=click.Path(path_type=Path), help="MetaTrawl cache directory containing genomes/ and genes/.")
@click.option("--output-dir", type=click.Path(path_type=Path), help="Optional legacy directory for unified matrix reference files.")
@click.option("--genome", help="Build files for one genome accession.")
@click.option("--accessions", "accessions_file", type=click.Path(path_type=Path), help="Optional accession list; defaults to every cached genome.")
def cache_sync_matrix_files(
    cache_dir: Path,
    output_dir: Path,
    genome: str | None,
    accessions_file: Path | None,
) -> None:
    """Refresh per-genome BED, STB, and gene-range files from the current cache."""
    if genome is not None and accessions_file is not None:
        raise click.UsageError("Use either --genome or --accessions, not both.")
    accessions = cache.read_accessions_file(accessions_file) if accessions_file is not None else None
    cache_dir = Path(cache_dir)
    try:
        synced = cache.sync_matrix_requirement_files(
            cache_dir=cache_dir,
            accessions=accessions,
            genome=genome,
        )
        legacy_files = (
            cache.build_matrix_reference_files(
                genome_dir=cache_dir / "genomes",
                gene_dir=cache_dir / "genes",
                output_dir=output_dir,
                accessions=accessions,
                genome=genome,
            )
            if output_dir is not None
            else None
        )
    except (FileNotFoundError, ValueError) as exc:
        raise click.ClickException(str(exc)) from exc
    payload = {
        "accessions": len(synced.accessions),
        "bed_dir": str(synced.bed_dir),
        "stb_dir": str(synced.stb_dir),
        "gene_range_dir": str(synced.gene_range_dir),
    }
    if legacy_files is not None:
        payload["legacy_reference_fasta"] = str(legacy_files.reference_fasta)
        payload["legacy_gene_fasta"] = str(legacy_files.gene_fasta)
        payload["legacy_bed_file"] = str(legacy_files.bed_file)
        payload["legacy_stb_file"] = str(legacy_files.stb_file)
        payload["legacy_gene_range_table"] = str(legacy_files.gene_range_table)
    click.echo(json.dumps(payload))


@cache_group.command("serve")
@click.option("--cache-dir", required=True, type=click.Path(path_type=Path), help="Shared genome cache directory.")
@click.option("--host", default="127.0.0.1", show_default=True)
@click.option("--port", default=8765, show_default=True, type=int)
def cache_serve(cache_dir: Path, host: str, port: int) -> None:
    """Serve cache prepare requests over a small local HTTP API."""
    logger = WorkflowLogger()
    manager = cache.GenomeCache(cache_dir, logger=logger)

    class Handler(BaseHTTPRequestHandler):
        def do_POST(self) -> None:  # noqa: N802 - stdlib handler API
            if self.path != "/prepare":
                self.send_error(404)
                return
            length = int(self.headers.get("Content-Length", "0"))
            payload = json.loads(self.rfile.read(length) or b"{}")
            try:
                prepared = manager.prepare_reference(
                    accessions=[str(item) for item in payload["accessions"]],
                    output_dir=Path(payload["output_dir"]),
                    sample=payload.get("sample"),
                )
                body = json.dumps({"reference_fasta": str(prepared.reference_fasta), "gene_fasta": str(prepared.gene_fasta)}).encode()
                self.send_response(200)
                self.send_header("Content-Type", "application/json")
                self.send_header("Content-Length", str(len(body)))
                self.end_headers()
                self.wfile.write(body)
            except Exception as exc:
                body = json.dumps({"error": str(exc)}).encode()
                self.send_response(500)
                self.send_header("Content-Type", "application/json")
                self.send_header("Content-Length", str(len(body)))
                self.end_headers()
                self.wfile.write(body)

        def log_message(self, format: str, *args: object) -> None:
            return

    click.echo(f"serving cache_dir={cache_dir} host={host} port={port}")
    ThreadingHTTPServer((host, port), Handler).serve_forever()


@cli.command("profile-sra")
@click.option("--db", "db_file", required=True, type=click.Path(path_type=Path), help="MetaTrawl DuckDB registry.")
@click.option("--remaining-csv", required=True, type=click.Path(path_type=Path), help="CSV from profiles remaining.")
@click.option("--cache-dir", required=True, type=click.Path(path_type=Path), help="Shared genome cache directory.")
@click.option("--scratch-dir", required=True, type=click.Path(path_type=Path), help="Disposable worker scratch directory.")
@click.option("--sylph-db", type=click.Path(path_type=Path), help="Sylph database used to select genomes from reads.")
@click.option("--output-dir", type=click.Path(path_type=Path), help="Directory where Sylph/profile outputs are written.")
@click.option("--accessions-dir", type=click.Path(path_type=Path), help="Manual override directory with per-run accession files.")
@click.option("--threads", type=int, default=8, show_default=True)
def profile_sra(
    db_file: Path,
    remaining_csv: Path,
    cache_dir: Path,
    scratch_dir: Path,
    sylph_db: Path | None,
    output_dir: Path | None,
    accessions_dir: Path | None,
    threads: int,
) -> None:
    """Run the local SRA profiling lifecycle with scratch cleanup."""
    run_ids = _read_remaining_csv(remaining_csv)
    workflows.profile_sra_runs(
        run_ids=run_ids,
        db_file=db_file,
        cache_dir=cache_dir,
        scratch_dir=scratch_dir,
        sylph_db=sylph_db,
        output_dir=output_dir,
        accessions_dir=accessions_dir,
        threads=threads,
        logger=WorkflowLogger(),
    )


@cli.command("sync-profile")
@click.option("--db", "db_file", required=True, type=click.Path(path_type=Path), help="MetaTrawl DuckDB registry.")
@click.option("--cache-dir", required=True, type=click.Path(path_type=Path), help="Shared genome cache directory.")
@click.option("--scratch-dir", required=True, type=click.Path(path_type=Path), help="Disposable worker scratch directory.")
@click.option("--output-dir", required=True, type=click.Path(path_type=Path), help="Directory where profile outputs are written.")
@click.option("--sylph-db", type=click.Path(path_type=Path), help="Sylph database used to select genomes from reads.")
@click.option("--accessions-dir", type=click.Path(path_type=Path), help="Manual override directory with per-run accession files.")
@click.option("--threads", type=int, default=8, show_default=True)
@click.option("--skip-dependency-check", is_flag=True, help="Do not preflight external tools before syncing.")
@click.option("--keep-profile-outputs", is_flag=True, help="Keep imported profile/stat/Sylph files on disk for debugging.")
def sync_profile(
    db_file: Path,
    cache_dir: Path,
    scratch_dir: Path,
    output_dir: Path,
    sylph_db: Path | None,
    accessions_dir: Path | None,
    threads: int,
    skip_dependency_check: bool,
    keep_profile_outputs: bool,
) -> None:
    """Profile all remaining runs and import completed outputs into DuckDB."""
    try:
        summary = workflows.sync_remaining_profiles(
            db_file=db_file,
            cache_dir=cache_dir,
            scratch_dir=scratch_dir,
            output_dir=output_dir,
            sylph_db=sylph_db,
            accessions_dir=accessions_dir,
            threads=threads,
            check_dependencies=not skip_dependency_check,
            cleanup_outputs=not keep_profile_outputs,
            logger=WorkflowLogger(),
        )
    except (FileNotFoundError, RuntimeError, ValueError, healthcheck.DependencyCheckError) as exc:
        raise click.ClickException(str(exc)) from exc
    click.echo(
        f"sync-profile requested={summary.requested} imported={summary.imported} "
        f"skipped={summary.skipped} failed={summary.failed} "
        f"cleaned_files={summary.cleaned_files}"
    )
    if summary.failed or summary.skipped:
        raise click.ClickException(
            "Sync completed with failures. Successful samples were checkpointed; "
            "rerun the same command to retry remaining samples."
        )


@cli.group("matrix")
def matrix_group() -> None:
    """Build and compare ZipStrain matrix stores from imported samples."""


@matrix_group.command("build")
@click.option("--db", "db_file", required=True, type=click.Path(path_type=Path), help="MetaTrawl DuckDB registry.")
@click.option("--genome", required=True, help="Genome scope for the matrix store.")
@click.option("--bed-file", required=True, type=click.Path(path_type=Path), help="BED file defining scaffold extents.")
@click.option("--stb-file", required=True, type=click.Path(path_type=Path), help="STB file defining scaffold-to-genome mapping.")
@click.option("--gene-range-table", type=click.Path(path_type=Path), help="Optional gene range table for gene ANI.")
@click.option("--output-file", required=True, type=click.Path(path_type=Path), help="Output matrix store.")
@click.option("--matrix-id", default=None, help="Registry ID for this matrix. Defaults to output-file stem.")
@click.option("--overwrite", is_flag=True, help="Replace an existing matrix file/registry row.")
@click.option("--sparse", is_flag=True, help="Store matrices using ZipStrain's sparse HDF5 matrix layout.")
@click.option("--min-coverage", type=float, default=None)
@click.option("--min-breadth", type=float, default=None)
@click.option("--min-ber", type=float, default=None)
@click.option("--min-sylph-abundance", type=float, default=None)
@click.option("--memory-limit-gb", type=float, default=16.0, show_default=True)
@click.option("--export-batch-mb", type=float, default=128.0, show_default=True)
def matrix_build(
    db_file: Path,
    genome: str,
    bed_file: Path,
    stb_file: Path,
    gene_range_table: Path | None,
    output_file: Path,
    matrix_id: str | None,
    overwrite: bool,
    sparse: bool,
    min_coverage: float | None,
    min_breadth: float | None,
    min_ber: float | None,
    min_sylph_abundance: float | None,
    memory_limit_gb: float,
    export_batch_mb: float,
) -> None:
    """Build a matrix store from filtered DuckDB profile rows."""
    filters = registry.MatrixFilters(
        min_coverage=min_coverage,
        min_breadth=min_breadth,
        min_ber=min_ber,
        min_sylph_abundance=min_sylph_abundance,
    )
    try:
        with registry.connect(db_file) as conn:
            store = workflows.build_matrix_from_database(
                conn,
                output_file=output_file,
                genome=genome,
                bed_file=bed_file,
                stb_file=stb_file,
                gene_range_table=gene_range_table,
                filters=filters,
                overwrite=overwrite,
                matrix_id=matrix_id,
                sparse=sparse,
                memory_limit_gb=memory_limit_gb,
                export_batch_mb=export_batch_mb,
                logger=WorkflowLogger(),
            )
    except (FileExistsError, FileNotFoundError, ValueError) as exc:
        raise click.ClickException(str(exc)) from exc
    click.echo(f"matrix_id={store.matrix_id} wrote={store.matrix_file} profiles={store.profile_count}")


@matrix_group.command("sync-build")
@click.option("--db", "db_file", required=True, type=click.Path(path_type=Path), help="MetaTrawl DuckDB registry.")
@click.option("--matrix-dir", required=True, type=click.Path(path_type=Path), help="Directory for per-genome HDF5 matrices.")
@click.option("--genome", "genomes", multiple=True, help="Optional genome accession/name to sync. Repeat for multiple genomes.")
@click.option("--bed-file", type=click.Path(path_type=Path), help="Optional cache-wide BED fallback file.")
@click.option("--stb-file", type=click.Path(path_type=Path), help="Optional cache-wide STB fallback file.")
@click.option("--bed-dir", type=click.Path(path_type=Path), help="Optional directory with per-genome ACCESSION.bed files.")
@click.option("--stb-dir", type=click.Path(path_type=Path), help="Optional directory with per-genome ACCESSION.stb files.")
@click.option("--gene-range-table", type=click.Path(path_type=Path), help="Optional cache-wide gene range table for gene ANI.")
@click.option("--gene-range-dir", type=click.Path(path_type=Path), help="Optional directory with per-genome ACCESSION.gene_ranges.tsv files.")
@click.option("--sparse", is_flag=True, help="Build missing matrices using ZipStrain's sparse HDF5 layout.")
@click.option("--min-coverage", type=float, default=None)
@click.option("--min-breadth", type=float, default=None)
@click.option("--min-ber", type=float, default=None)
@click.option("--min-sylph-abundance", type=float, default=None)
@click.option("--memory-limit-gb", type=float, default=16.0, show_default=True)
@click.option("--export-batch-mb", type=float, default=128.0, show_default=True)
def matrix_sync_build(
    db_file: Path,
    matrix_dir: Path,
    genomes: tuple[str, ...],
    bed_file: Path,
    stb_file: Path,
    bed_dir: Path | None,
    stb_dir: Path | None,
    gene_range_table: Path | None,
    gene_range_dir: Path | None,
    sparse: bool,
    min_coverage: float | None,
    min_breadth: float | None,
    min_ber: float | None,
    min_sylph_abundance: float | None,
    memory_limit_gb: float,
    export_batch_mb: float,
) -> None:
    """Build missing per-genome matrices and append newly imported samples."""
    filters = registry.MatrixFilters(
        min_coverage=min_coverage,
        min_breadth=min_breadth,
        min_ber=min_ber,
        min_sylph_abundance=min_sylph_abundance,
    )
    try:
        with registry.connect(db_file) as conn:
            summary = workflows.sync_build_matrices(
                conn,
                matrix_dir=matrix_dir,
                genomes=list(genomes) if genomes else None,
                bed_file=bed_file,
                stb_file=stb_file,
                bed_dir=bed_dir,
                stb_dir=stb_dir,
                gene_range_table=gene_range_table,
                gene_range_dir=gene_range_dir,
                filters=filters,
                sparse=sparse,
                memory_limit_gb=memory_limit_gb,
                export_batch_mb=export_batch_mb,
                logger=WorkflowLogger(),
            )
    except (FileNotFoundError, RuntimeError, ValueError) as exc:
        raise click.ClickException(str(exc)) from exc
    click.echo(
        f"matrix-sync-build genomes={summary.genomes} built={summary.built} "
        f"appended={summary.appended} up_to_date={summary.up_to_date} "
        f"skipped={summary.skipped} failed={summary.failed}"
    )
    if summary.failed:
        raise click.ClickException("Matrix sync build completed with failures. Check the genome-level logs above.")


@matrix_group.command("append")
@click.option("--db", "db_file", required=True, type=click.Path(path_type=Path), help="MetaTrawl DuckDB registry.")
@click.option("--matrix-id", help="Registered matrix ID.")
@click.option("--matrix-file", type=click.Path(path_type=Path), help="Registered matrix HDF5 file.")
@click.option("--memory-limit-gb", type=float, default=16.0, show_default=True)
@click.option("--export-batch-mb", type=float, default=128.0, show_default=True)
def matrix_append(db_file: Path, matrix_id: str | None, matrix_file: Path | None, memory_limit_gb: float, export_batch_mb: float) -> None:
    """Append newly imported complete samples to a registered matrix store."""
    try:
        with registry.connect(db_file) as conn:
            resolved_matrix_file, resolved_matrix_id = workflows.resolve_matrix_file(
                conn,
                matrix_id=matrix_id,
                matrix_file=matrix_file,
            )
            appended = workflows.append_matrix_from_database(
                conn,
                matrix_file=resolved_matrix_file,
                matrix_id=resolved_matrix_id,
                memory_limit_gb=memory_limit_gb,
                export_batch_mb=export_batch_mb,
                logger=WorkflowLogger(),
            )
    except (FileNotFoundError, ValueError) as exc:
        raise click.ClickException(str(exc)) from exc
    click.echo(f"matrix_file={resolved_matrix_file} appended={appended}")


@matrix_group.command("compare")
@click.option("--db", "db_file", required=True, type=click.Path(path_type=Path), help="MetaTrawl DuckDB registry.")
@click.option("--matrix-id", help="Registered matrix ID.")
@click.option("--matrix-file", type=click.Path(path_type=Path), help="Registered matrix HDF5 file.")
@click.option("--output-file", required=True, type=click.Path(path_type=Path), help="Output ZipStrain compare DuckDB.")
@click.option("--calculate", default="all", show_default=True, help="Metrics to calculate.")
@click.option("--genome", default="all", show_default=True, help="Optional compare genome scope.")
@click.option("--backend", default="numpy", show_default=True, help="ZipStrain matrix backend.")
@click.option("--memory-limit-gb", type=float, default=16.0, show_default=True)
def matrix_compare(db_file: Path, matrix_id: str | None, matrix_file: Path | None, output_file: Path, calculate: str, genome: str, backend: str, memory_limit_gb: float) -> None:
    """Compare a registered matrix store."""
    try:
        with registry.connect(db_file) as conn:
            resolved_matrix_file, resolved_matrix_id = workflows.resolve_matrix_file(
                conn,
                matrix_id=matrix_id,
                matrix_file=matrix_file,
            )
            compare_id = workflows.compare_matrix_file(
                conn,
                matrix_file=resolved_matrix_file,
                matrix_id=resolved_matrix_id,
                output_file=output_file,
                calculate=calculate,
                genome=genome,
                backend=backend,
                memory_limit_gb=memory_limit_gb,
                logger=WorkflowLogger(),
            )
    except (FileNotFoundError, ValueError) as exc:
        raise click.ClickException(str(exc)) from exc
    click.echo(f"compare_id={compare_id} matrix_file={resolved_matrix_file} wrote={output_file}")


@matrix_group.command("sync-compare")
@click.option("--db", "db_file", required=True, type=click.Path(path_type=Path), help="MetaTrawl DuckDB registry.")
@click.option("--matrix-dir", required=True, type=click.Path(path_type=Path), help="Directory containing per-genome HDF5 matrices.")
@click.option("--compare-dir", required=True, type=click.Path(path_type=Path), help="Directory for per-genome compare DuckDB files.")
@click.option("--calculate", default="all", show_default=True, help="Metrics to calculate.")
@click.option("--genome", default="all", show_default=True, help="Optional compare genome scope.")
@click.option("--backend", default="numpy", show_default=True, help="ZipStrain matrix backend.")
@click.option("--memory-limit-gb", type=float, default=16.0, show_default=True)
def matrix_sync_compare(
    db_file: Path,
    matrix_dir: Path,
    compare_dir: Path,
    calculate: str,
    genome: str,
    backend: str,
    memory_limit_gb: float,
) -> None:
    """Run resumable compare for every matrix file in a directory."""
    try:
        with registry.connect(db_file) as conn:
            summary = workflows.sync_compare_matrices(
                conn,
                matrix_dir=matrix_dir,
                compare_dir=compare_dir,
                calculate=calculate,
                genome=genome,
                backend=backend,
                memory_limit_gb=memory_limit_gb,
                logger=WorkflowLogger(),
            )
    except (FileNotFoundError, RuntimeError, ValueError) as exc:
        raise click.ClickException(str(exc)) from exc
    click.echo(
        f"matrix-sync-compare matrices={summary.matrices} "
        f"compared={summary.compared} failed={summary.failed}"
    )
    if summary.failed:
        raise click.ClickException("Matrix sync compare completed with failures. Check the matrix-level logs above.")


@cli.command("status")
@click.option("--db", "db_file", required=True, type=click.Path(path_type=Path), help="MetaTrawl DuckDB registry.")
def status(db_file: Path) -> None:
    """Show high-level registry counts."""
    with registry.connect(db_file) as conn:
        counts = registry.registry_status(conn)
    for key, value in counts.items():
        click.echo(f"{key}={value}")


@cli.command("test")
def test() -> None:
    """Report installed tools and Python packages used by MetaTrawl workflows."""
    healthcheck.render_checks(healthcheck.collect_checks(), console=Console())


@cli.command("check")
def check() -> None:
    """Fail unless all dependencies needed by MetaTrawl are available."""
    checks = healthcheck.collect_checks()
    healthcheck.render_checks(checks, console=Console())
    try:
        healthcheck.assert_dependencies(checks)
    except healthcheck.DependencyCheckError as exc:
        raise click.ClickException(str(exc)) from exc


def _write_remaining_csv(output_file: Path, run_ids: list[str]) -> None:
    output_file.parent.mkdir(parents=True, exist_ok=True)
    with output_file.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["run_id"])
        writer.writeheader()
        for run_id in run_ids:
            writer.writerow({"run_id": run_id})


def _read_remaining_csv(input_file: Path) -> list[str]:
    with input_file.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None or "run_id" not in reader.fieldnames:
            raise click.UsageError("remaining CSV must contain a run_id column")
        return [row["run_id"] for row in reader if row.get("run_id")]


def _read_profile_manifest(manifest: Path) -> list[registry.ProfileBundle]:
    required = {"run_id", "profile_file", "genome_stats_file", "sylph_abundance_file"}
    with manifest.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise click.UsageError("Profile manifest is empty.")
        missing = required - set(reader.fieldnames)
        if missing:
            raise click.UsageError(f"Profile manifest missing required columns: {', '.join(sorted(missing))}")
        return [
            registry.ProfileBundle(
                run_id=row["run_id"],
                profile_file=Path(row["profile_file"]),
                genome_stats_file=Path(row["genome_stats_file"]),
                gene_stats_file=Path(row["gene_stats_file"]) if row.get("gene_stats_file") else None,
                sylph_abundance_file=Path(row["sylph_abundance_file"]),
            )
            for row in reader
        ]


def _print_rows(rows: list[dict[str, object]], *, columns: list[str], title: str) -> None:
    console = Console()
    table = Table(title=title)
    for column in columns:
        table.add_column(column)
    for row in rows:
        table.add_row(*(str(row.get(column, "")) for column in columns))
    console.print(table)


if __name__ == "__main__":
    cli()

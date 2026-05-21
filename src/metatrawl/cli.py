"""MetaTrawl command-line interface."""

from __future__ import annotations

import csv
from pathlib import Path

import click
from rich.console import Console
from rich.table import Table

from metatrawl import __version__
from metatrawl import db as registry
from metatrawl import healthcheck
from metatrawl import workflows


@click.group()
@click.version_option(version=__version__, prog_name="metatrawl")
def cli() -> None:
    """Registry tooling for SRA-scale ZipStrain profiling projects."""


@cli.command("init")
@click.option("--db", "db_file", required=True, type=click.Path(path_type=Path), help="MetaTrawl DuckDB registry.")
def init(db_file: Path) -> None:
    """Initialize a MetaTrawl registry."""
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
    _print_rows(rows, columns=["run_id", "deleted_at"], title="SRA runs")


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
    """Manage completed ZipStrain/Sylph profile bundles."""


@profiles_group.command("remaining")
@click.option("--db", "db_file", required=True, type=click.Path(path_type=Path), help="MetaTrawl DuckDB registry.")
@click.option("--output-file", "-o", type=click.Path(path_type=Path), help="CSV path for remaining SRA run IDs.")
def profiles_remaining(db_file: Path, output_file: Path | None) -> None:
    """Write active SRA run IDs that still need profile bundles."""
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
        columns=["run_id", "profile_file", "genome_stats_file", "sylph_abundance_file"],
        title="Profiles",
    )


@cli.group("matrix")
def matrix_group() -> None:
    """Build and compare ZipStrain matrix stores from registered profiles."""


@matrix_group.command("build")
@click.option("--db", "db_file", required=True, type=click.Path(path_type=Path), help="MetaTrawl DuckDB registry.")
@click.option("--genome", required=True, help="Genome scope for the matrix store.")
@click.option("--output-file", required=True, type=click.Path(path_type=Path), help="Output matrix store.")
@click.option("--matrix-id", default=None, help="Registry ID for this matrix. Defaults to output-file stem.")
@click.option("--overwrite", is_flag=True, help="Replace an existing matrix file/registry row.")
@click.option("--memory-limit-gb", type=float, default=16.0, show_default=True)
@click.option("--export-batch-mb", type=float, default=128.0, show_default=True)
def matrix_build(
    db_file: Path,
    genome: str,
    output_file: Path,
    matrix_id: str | None,
    overwrite: bool,
    memory_limit_gb: float,
    export_batch_mb: float,
) -> None:
    """Build a matrix store from complete registered profiles."""
    try:
        with registry.connect(db_file) as conn:
            store = workflows.build_matrix_from_registry(
                conn,
                output_file=output_file,
                genome=genome,
                overwrite=overwrite,
                matrix_id=matrix_id,
                memory_limit_gb=memory_limit_gb,
                export_batch_mb=export_batch_mb,
            )
    except (FileExistsError, ValueError) as exc:
        raise click.ClickException(str(exc)) from exc
    click.echo(f"matrix_id={store.matrix_id} wrote={store.matrix_file} profiles={store.profile_count}")


@matrix_group.command("compare")
@click.option("--db", "db_file", required=True, type=click.Path(path_type=Path), help="MetaTrawl DuckDB registry.")
@click.option("--matrix-id", required=True, help="Registered matrix ID.")
@click.option("--output-file", required=True, type=click.Path(path_type=Path), help="Output ZipStrain compare DuckDB.")
@click.option("--calculate", default="all", show_default=True, help="Metrics to calculate.")
@click.option("--genome", default="all", show_default=True, help="Optional compare genome scope.")
@click.option("--backend", default="numpy", show_default=True, help="ZipStrain matrix backend.")
@click.option("--memory-limit-gb", type=float, default=16.0, show_default=True)
def matrix_compare(
    db_file: Path,
    matrix_id: str,
    output_file: Path,
    calculate: str,
    genome: str,
    backend: str,
    memory_limit_gb: float,
) -> None:
    """Compare a registered matrix store."""
    try:
        with registry.connect(db_file) as conn:
            compare_id = workflows.compare_matrix_from_registry(
                conn,
                matrix_id=matrix_id,
                output_file=output_file,
                calculate=calculate,
                genome=genome,
                backend=backend,
                memory_limit_gb=memory_limit_gb,
            )
    except ValueError as exc:
        raise click.ClickException(str(exc)) from exc
    click.echo(f"compare_id={compare_id} wrote={output_file}")


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


def _write_remaining_csv(output_file: Path, run_ids: list[str]) -> None:
    output_file.parent.mkdir(parents=True, exist_ok=True)
    with output_file.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["run_id"])
        writer.writeheader()
        for run_id in run_ids:
            writer.writerow({"run_id": run_id})


def _read_profile_manifest(manifest: Path) -> list[registry.ProfileBundle]:
    required = {"run_id", "profile_file", "genome_stats_file", "sylph_abundance_file"}
    with manifest.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise click.UsageError("Profile manifest is empty.")
        missing = required - set(reader.fieldnames)
        if missing:
            raise click.UsageError(f"Profile manifest missing required columns: {', '.join(sorted(missing))}")
        bundles = [
            registry.ProfileBundle(
                run_id=row["run_id"],
                profile_file=Path(row["profile_file"]),
                genome_stats_file=Path(row["genome_stats_file"]),
                sylph_abundance_file=Path(row["sylph_abundance_file"]),
            )
            for row in reader
        ]
    return bundles


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

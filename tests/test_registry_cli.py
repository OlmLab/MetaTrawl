from pathlib import Path
import csv
import types
import sys

from click.testing import CliRunner
import duckdb

from metatrawl import cli


def _touch(path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("x")
    return path


def _write_manifest(path: Path, rows: list[dict[str, str]], fieldnames: list[str] | None = None) -> Path:
    fieldnames = fieldnames or ["run_id", "profile_file", "genome_stats_file", "sylph_abundance_file"]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    return path


def _add_complete_profile(runner: CliRunner, db_file: Path, tmp_path: Path, run_id: str) -> None:
    profile = _touch(tmp_path / f"{run_id}.parquet")
    genome_stats = _touch(tmp_path / f"{run_id}_genome_stats.parquet")
    sylph = _touch(tmp_path / f"{run_id}_sylph.tsv")
    manifest = _write_manifest(
        tmp_path / f"{run_id}_manifest.csv",
        [
            {
                "run_id": run_id,
                "profile_file": str(profile),
                "genome_stats_file": str(genome_stats),
                "sylph_abundance_file": str(sylph),
            }
        ],
    )
    result = runner.invoke(cli.cli, ["profiles", "add", "--db", str(db_file), "--manifest", str(manifest)])
    assert result.exit_code == 0, result.output


def test_init_creates_duckdb_schema(tmp_path: Path):
    runner = CliRunner()
    db_file = tmp_path / "metatrawl.duckdb"

    result = runner.invoke(cli.cli, ["init", "--db", str(db_file)])

    assert result.exit_code == 0, result.output
    tables = {
        row[0]
        for row in duckdb.connect(str(db_file))
        .execute("SHOW TABLES")
        .fetchall()
    }
    assert {"sra_runs", "profiles", "matrix_stores", "matrix_compares"} <= tables


def test_adding_run_ids_is_idempotent(tmp_path: Path):
    runner = CliRunner()
    db_file = tmp_path / "metatrawl.duckdb"

    first = runner.invoke(cli.cli, ["runs", "add", "--db", str(db_file), "SRR1", "SRR1", "SRR2"])
    second = runner.invoke(cli.cli, ["runs", "add", "--db", str(db_file), "SRR1", "SRR2"])

    assert first.exit_code == 0, first.output
    assert "added=2" in first.output
    assert second.exit_code == 0, second.output
    assert "added=0" in second.output


def test_deleted_run_is_excluded_from_remaining_profiles(tmp_path: Path):
    runner = CliRunner()
    db_file = tmp_path / "metatrawl.duckdb"

    assert runner.invoke(cli.cli, ["runs", "add", "--db", str(db_file), "SRR1", "SRR2"]).exit_code == 0
    assert runner.invoke(cli.cli, ["runs", "delete", "--db", str(db_file), "SRR2"]).exit_code == 0
    output = tmp_path / "remaining.csv"
    result = runner.invoke(cli.cli, ["profiles", "remaining", "--db", str(db_file), "--output-file", str(output)])

    assert result.exit_code == 0, result.output
    assert output.read_text().splitlines() == ["run_id", "SRR1"]


def test_profiles_remaining_emits_only_unprofiled_runs(tmp_path: Path):
    runner = CliRunner()
    db_file = tmp_path / "metatrawl.duckdb"
    assert runner.invoke(cli.cli, ["runs", "add", "--db", str(db_file), "SRR1", "SRR2"]).exit_code == 0
    _add_complete_profile(runner, db_file, tmp_path, "SRR1")

    output = tmp_path / "remaining.csv"
    result = runner.invoke(cli.cli, ["profiles", "remaining", "--db", str(db_file), "--output-file", str(output)])

    assert result.exit_code == 0, result.output
    assert output.read_text().splitlines() == ["run_id", "SRR2"]


def test_profiles_remaining_reports_all_complete_and_writes_header(tmp_path: Path):
    runner = CliRunner()
    db_file = tmp_path / "metatrawl.duckdb"
    assert runner.invoke(cli.cli, ["runs", "add", "--db", str(db_file), "SRR1"]).exit_code == 0
    _add_complete_profile(runner, db_file, tmp_path, "SRR1")

    output = tmp_path / "remaining.csv"
    result = runner.invoke(cli.cli, ["profiles", "remaining", "--db", str(db_file), "--output-file", str(output)])

    assert result.exit_code == 0, result.output
    assert "All added runs have complete profiles" in result.output
    assert output.read_text().splitlines() == ["run_id"]


def test_profiles_add_imports_completed_profile_bundle(tmp_path: Path):
    runner = CliRunner()
    db_file = tmp_path / "metatrawl.duckdb"
    assert runner.invoke(cli.cli, ["runs", "add", "--db", str(db_file), "SRR1"]).exit_code == 0

    _add_complete_profile(runner, db_file, tmp_path, "SRR1")
    result = runner.invoke(cli.cli, ["profiles", "list", "--db", str(db_file)])

    assert result.exit_code == 0, result.output
    assert "SRR1" in result.output
    rows = duckdb.connect(str(db_file)).execute("SELECT profile_file FROM profiles WHERE run_id = 'SRR1'").fetchall()
    assert Path(rows[0][0]).name == "SRR1.parquet"


def test_profiles_add_rejects_missing_required_columns(tmp_path: Path):
    runner = CliRunner()
    db_file = tmp_path / "metatrawl.duckdb"
    manifest = _write_manifest(
        tmp_path / "bad_manifest.csv",
        [{"run_id": "SRR1", "profile_file": "profile.parquet"}],
        fieldnames=["run_id", "profile_file"],
    )

    result = runner.invoke(cli.cli, ["profiles", "add", "--db", str(db_file), "--manifest", str(manifest)])

    assert result.exit_code != 0
    assert "missing required columns" in result.output


def test_profiles_add_rejects_unknown_runs_unless_add_runs_is_used(tmp_path: Path):
    runner = CliRunner()
    db_file = tmp_path / "metatrawl.duckdb"
    profile = _touch(tmp_path / "SRR1.parquet")
    genome_stats = _touch(tmp_path / "SRR1_genome_stats.parquet")
    sylph = _touch(tmp_path / "SRR1_sylph.tsv")
    manifest = _write_manifest(
        tmp_path / "manifest.csv",
        [
            {
                "run_id": "SRR1",
                "profile_file": str(profile),
                "genome_stats_file": str(genome_stats),
                "sylph_abundance_file": str(sylph),
            }
        ],
    )

    rejected = runner.invoke(cli.cli, ["profiles", "add", "--db", str(db_file), "--manifest", str(manifest)])
    accepted = runner.invoke(
        cli.cli,
        ["profiles", "add", "--db", str(db_file), "--manifest", str(manifest), "--add-runs"],
    )

    assert rejected.exit_code != 0
    assert accepted.exit_code == 0, accepted.output
    assert "imported=1" in accepted.output


def test_matrix_build_uses_only_complete_profiles(tmp_path: Path, monkeypatch):
    runner = CliRunner()
    db_file = tmp_path / "metatrawl.duckdb"
    assert runner.invoke(cli.cli, ["runs", "add", "--db", str(db_file), "SRR1", "SRR2"]).exit_code == 0
    _add_complete_profile(runner, db_file, tmp_path, "SRR1")
    calls: list[dict[str, object]] = []

    def build_matrix_hdf5(**kwargs):
        calls.append({
            **kwargs,
            "staged_profile_count": len(list(Path(kwargs["profile_dir"]).glob("*.parquet"))),
        })
        Path(kwargs["output_file"]).write_text("matrix")

    zipstrain_module = types.ModuleType("zipstrain")
    matrix_pairs_module = types.ModuleType("zipstrain.matrix_pairs")
    matrix_pairs_module.build_matrix_hdf5 = build_matrix_hdf5
    monkeypatch.setitem(sys.modules, "zipstrain", zipstrain_module)
    monkeypatch.setitem(sys.modules, "zipstrain.matrix_pairs", matrix_pairs_module)
    monkeypatch.setattr(zipstrain_module, "matrix_pairs", matrix_pairs_module, raising=False)

    result = runner.invoke(
        cli.cli,
        [
            "matrix",
            "build",
            "--db",
            str(db_file),
            "--genome",
            "genome_a",
            "--output-file",
            str(tmp_path / "matrix.h5"),
        ],
    )

    assert result.exit_code == 0, result.output
    assert "profiles=1" in result.output
    assert len(calls) == 1
    assert calls[0]["staged_profile_count"] == 1


def test_matrix_build_fails_when_no_complete_profiles_exist(tmp_path: Path):
    runner = CliRunner()
    db_file = tmp_path / "metatrawl.duckdb"
    assert runner.invoke(cli.cli, ["runs", "add", "--db", str(db_file), "SRR1"]).exit_code == 0

    result = runner.invoke(
        cli.cli,
        [
            "matrix",
            "build",
            "--db",
            str(db_file),
            "--genome",
            "genome_a",
            "--output-file",
            str(tmp_path / "matrix.h5"),
        ],
    )

    assert result.exit_code != 0
    assert "No complete profiles" in result.output


def test_matrix_compare_uses_registered_matrix_store(tmp_path: Path, monkeypatch):
    runner = CliRunner()
    db_file = tmp_path / "metatrawl.duckdb"
    matrix_file = _touch(tmp_path / "matrix.h5")
    calls: list[dict[str, object]] = []

    with duckdb.connect(str(db_file)) as conn:
        from metatrawl import db as registry

        registry.init_schema(conn)
        registry.register_matrix_store(
            conn,
            matrix_id="genome_a",
            genome="genome_a",
            matrix_file=matrix_file,
            profile_count=2,
        )

    def matrix_compare(**kwargs):
        calls.append(kwargs)
        Path(kwargs["output_file"]).write_text("compare")

    zipstrain_module = types.ModuleType("zipstrain")
    matrix_pairs_module = types.ModuleType("zipstrain.matrix_pairs")
    matrix_pairs_module.matrix_compare = matrix_compare
    monkeypatch.setitem(sys.modules, "zipstrain", zipstrain_module)
    monkeypatch.setitem(sys.modules, "zipstrain.matrix_pairs", matrix_pairs_module)
    monkeypatch.setattr(zipstrain_module, "matrix_pairs", matrix_pairs_module, raising=False)

    result = runner.invoke(
        cli.cli,
        [
            "matrix",
            "compare",
            "--db",
            str(db_file),
            "--matrix-id",
            "genome_a",
            "--output-file",
            str(tmp_path / "compare.duckdb"),
            "--calculate",
            "all",
        ],
    )

    assert result.exit_code == 0, result.output
    assert "compare_id=compare" in result.output
    assert len(calls) == 1
    assert calls[0]["matrix_db_file"] == matrix_file
    with duckdb.connect(str(db_file)) as conn:
        rows = conn.execute("SELECT compare_id, matrix_id, calculate FROM matrix_compares").fetchall()
    assert rows == [("compare", "genome_a", "all")]

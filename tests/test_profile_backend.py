"""Profiling backend selection belongs to execution, not scientific filters."""

from dataclasses import replace
import json
from pathlib import Path
import shlex
import subprocess

from click.testing import CliRunner
import pytest

from metatrawl import cache, cli, db, provenance, workflows
from metatrawl.config import ProfileConfig, WorkflowConfig, load_workflow_config
from metatrawl.execution import WorkflowRuntime
from metatrawl.logging import WorkflowLogger
from test_registry_cli import _write_bundle_files


@pytest.mark.parametrize("suffix", ["toml", "json"])
@pytest.mark.parametrize("backend", ["python", "rust_profiler", "custom_plugin"])
def test_backend_config_accepts_zipstrain_plugin_names(tmp_path, suffix, backend):
    config_file = tmp_path / f"workflow.{suffix}"
    config_file.write_text(
        f'[profile]\nbackend = "{backend}"\n' if suffix == "toml"
        else json.dumps({"profile": {"backend": backend}})
    )
    config = load_workflow_config(config_file, threads=4, sample_count=2)
    assert config.profile.backend == backend
    assert config.profile.min_freq == 0.05
    assert config.matrix_compare.backend is None


def test_default_backend_remains_python():
    assert load_workflow_config(None, threads=4, sample_count=2).profile.backend == "python"


@pytest.mark.parametrize("value", ["", "   ", "rust profiler", 3, False, None])
def test_invalid_backend_rejected_in_config(tmp_path, value):
    path = tmp_path / "workflow.json"
    path.write_text(json.dumps({"profile": {"backend": value}}))
    with pytest.raises(ValueError, match="profile.backend must be a non-empty backend name"):
        load_workflow_config(path, threads=1, sample_count=1)


@pytest.mark.parametrize("command", ["sync-profile", "profile-sra"])
@pytest.mark.parametrize("cli_backend,toml_backend,expected", [
    (None, None, "python"), (None, "rust_profiler", "rust_profiler"),
    ("rust_profiler", "python", "rust_profiler"), ("python", "rust_profiler", "python"),
])
def test_cli_backend_precedence(tmp_path, monkeypatch, command, cli_backend, toml_backend, expected):
    selected = []

    def sync(**kwargs):
        selected.append(kwargs["workflow_config"])
        return workflows.SyncSummary(requested=0, imported=0, skipped=0, failed=0, cleaned_files=0)

    def profile(**kwargs):
        selected.append(kwargs["workflow_config"])
        return {}

    monkeypatch.setattr(workflows, "sync_remaining_profiles", sync)
    monkeypatch.setattr(workflows, "profile_sra_runs", profile)
    args = [
        command, "--db", str(tmp_path / "project.duckdb"),
        "--cache-dir", str(tmp_path / "cache"), "--scratch-dir", str(tmp_path / "scratch"),
        "--output-dir", str(tmp_path / "outputs"),
    ]
    if command == "profile-sra":
        remaining = tmp_path / "remaining.csv"
        remaining.write_text("run_id\nSRR1\n")
        args.extend(["--remaining-csv", str(remaining)])
    if toml_backend is not None:
        config = tmp_path / "workflow.toml"
        config.write_text(f'[profile]\nbackend = "{toml_backend}"\n[stages.profile]\nexecution = "slurm"\nthreads = 12\n')
        args.extend(["--workflow-config", str(config)])
    if cli_backend is not None:
        args.extend(["--backend", cli_backend])
    result = CliRunner().invoke(cli.cli, args)
    assert result.exit_code == 0, result.output
    assert selected[0].profile.backend == expected
    if toml_backend is not None:
        assert selected[0].stage("profile").execution == "slurm"
        assert selected[0].stage("profile").threads == 12


def test_empty_cli_backend_fails_cleanly(tmp_path):
    result = CliRunner().invoke(cli.cli, [
        "sync-profile", "--db", str(tmp_path / "db.duckdb"),
        "--cache-dir", "cache", "--scratch-dir", "scratch", "--output-dir", "outputs",
        "--backend", "",
    ])
    assert result.exit_code != 0
    assert "profile.backend must be a non-empty backend name" in result.output


def test_backend_switch_preserves_contract_and_accepts_import(tmp_path):
    python = ProfileConfig()
    rust = replace(python, backend="rust_profiler")
    assert provenance.profiling_contract(python) == provenance.profiling_contract(rust)
    config = replace(WorkflowConfig.legacy(threads=2, sample_count=1), profile=rust)
    assert provenance.operational_config(config)["profile"]["backend"] == "rust_profiler"
    with db.connect(tmp_path / "project.duckdb") as conn:
        db.add_runs(conn, ["PYTHON", "RUST"])
        db.import_profile_bundle(conn, _write_bundle_files(tmp_path, "PYTHON"))
        bundle = _write_bundle_files(tmp_path, "RUST")
        provenance.write_manifest(
            profile_file=bundle.profile_file, sample_id="RUST",
            contract=provenance.profiling_contract(rust), references={},
            details={"profile_backend": "rust_profiler"},
        )
        db.import_profile_bundle(conn, bundle)
        assert conn.execute("SELECT count(*) FROM profiling_contracts").fetchone()[0] == 1
        assert conn.execute("SELECT count(*) FROM samples WHERE status='complete'").fetchone()[0] == 2


@pytest.mark.parametrize("backend", ["python", "rust_profiler"])
@pytest.mark.parametrize("execution", ["local", "slurm"])
def test_profile_backend_reaches_worker_command_and_manifest(tmp_path, monkeypatch, backend, execution, capsys):
    sample = tmp_path / "scratch" / "SRR1"
    work = sample / "zipstrain_profile"
    work.mkdir(parents=True)
    config = replace(ProfileConfig(), backend=backend, min_freq=0.02)
    # The same scientific checkpoint was made using Python before this run.
    previous_config = replace(config, backend="python")
    (sample / "profiling_settings.json").write_text(provenance.canonical(provenance.profiling_contract(previous_config)))
    (sample / "SRR1.bam").write_bytes(b"test-bam")
    (sample / "read_layout.txt").write_text("paired\n")
    reference_dir = tmp_path / "reference"
    reference_dir.mkdir()
    fasta, genes, stb = reference_dir / "reference.fna", reference_dir / "genes.fna", reference_dir / "reference.stb"
    fasta.write_text(">chr\nACGT\n")
    genes.write_text(">chr_1 # 1 # 4 # 1 # ID=1\nACGT\n")
    stb.write_text("chr\tgenome\n")
    for name in ("genomes_bed_file.bed", "gene_range_table.tsv", "null_model.parquet", "reference.fasta"):
        (work / name).write_text("test-fixture")
    (work / "profiling_contract.json").write_text("{}")
    monkeypatch.setattr(workflows, "_validate_bam", lambda *args, **kwargs: None)
    monkeypatch.setattr(workflows, "_ensure_reference_fai", lambda **kwargs: None)
    captured = []
    scripts = []

    def write_outputs(command):
        captured.append(command)
        for suffix in ("profile", "genome_stats", "gene_stats"):
            (work / f"SRR1_{suffix}.parquet").write_text("fixture")

    def runner(command, **kwargs):
        if command[0] == "sbatch":
            script = Path(command[-1]).read_text()
            scripts.append(script)
            # Inspect the real generated job, then simulate successful outputs.
            write_outputs(shlex.split(script.splitlines()[-1]))
        else:
            write_outputs(list(command))
        return subprocess.CompletedProcess(command, 0, stdout="123\n", stderr="")

    workflow = replace(WorkflowConfig.legacy(threads=3, sample_count=1), profile=config)
    workflow = replace(workflow, stages={**workflow.stages, "profile": replace(workflow.stage("profile"), execution=execution)})
    runtime = WorkflowRuntime(workflow, state_dir=tmp_path / "state", logger=WorkflowLogger(), runner=runner)
    outputs = tmp_path / "outputs"
    outputs.mkdir()
    (outputs / "SRR1.sylph.tsv").write_text("Genome_file\tabundance\nGCF_1\t1\n")
    workflows._run_alignment_and_profile(
        run_id="SRR1", sample_scratch=sample,
        reference=cache.PreparedReference(reference_fasta=fasta, gene_fasta=genes, stb_file=stb),
        output_dir=outputs, runtime=runtime, profile_config=config, logger=WorkflowLogger(),
    )
    command = captured[0]
    assert command[:3] == ["zipstrain", "utilities", "profile-single"]
    assert command[command.index("--min-freq") + 1] == "0.02"
    assert command[command.index("--max-concurrency") + 1] == "3"
    if backend == "python":
        assert "--backend" not in command
    else:
        assert command[command.index("--backend") + 1] == backend
    assert bool(scripts) == (execution == "slurm")
    manifest = json.loads(provenance.manifest_path(outputs / "SRR1.profile.parquet").read_text())
    assert manifest["details"]["profile_backend"] == backend
    assert "backend" not in manifest["contract"]["profile"]
    assert f"backend={backend}" in capsys.readouterr().err

"""Canonical and legacy kernel selections must work across local/Slurm routes."""

from pathlib import Path
import json
import shlex
import subprocess
from types import SimpleNamespace

from click.testing import CliRunner
import duckdb
import polars as pl
import pytest

from metatrawl import cli, config, db, workflows
from metatrawl.execution import WorkflowRuntime
from metatrawl.logging import WorkflowLogger
from zipstrain import matrix_pairs as mp


@pytest.mark.parametrize("extension", ["toml", "json"])
@pytest.mark.parametrize("kernel,expected", [
    ("anchor", "anchor"), ("block", "block"), ("pairwise", "anchor"), ("gemm", "block"),
])
def test_kernel_configuration(tmp_path, extension, kernel, expected):
    path = tmp_path / f"workflow.{extension}"
    path.write_text(
        f'[matrix_compare]\ncompare_kernel="{kernel}"\n'
        if extension == "toml" else json.dumps({"matrix_compare": {"compare_kernel": kernel}})
    )
    parsed = config.load_workflow_config(path, threads=1, sample_count=1)
    assert parsed.matrix_compare.compare_kernel == expected
    assert parsed.matrix_compare.kwargs() == {"compare_kernel": expected}


def test_kernel_is_optional_and_invalid_values_are_rejected(tmp_path):
    assert config.load_workflow_config(None, threads=1, sample_count=1).matrix_compare.compare_kernel is None
    with pytest.raises(ValueError, match="compare_kernel"):
        config.MatrixCompareConfig(compare_kernel="invalid")
    path = tmp_path / "workflow.toml"
    path.write_text('[matrix_compare]\ncompare_kernel="invalid"\n')
    with pytest.raises(ValueError, match="compare_kernel"):
        config.load_workflow_config(path, threads=1, sample_count=1)


@pytest.mark.parametrize("command", ["compare", "sync-compare"])
@pytest.mark.parametrize("configured,override,expected", [
    (None, None, "anchor"), ("gemm", None, "block"),
    ("block", "pairwise", "anchor"), ("pairwise", "block", "block"),
])
def test_cli_and_toml_kernel_precedence(tmp_path, monkeypatch, command, configured, override, expected):
    matrices = tmp_path / "matrices"; matrices.mkdir()
    matrix = matrices / "genome.h5"; matrix.write_text("fixture")
    settings = tmp_path / "workflow.toml"
    settings.write_text('[matrix_compare]\ncalculate="genome_ani"\n' +
                        (f'compare_kernel="{configured}"\n' if configured else ""))
    calls = []
    def compare(conn, **kwargs):
        calls.append(kwargs)
        return "genome"
    monkeypatch.setattr(workflows, "compare_matrix_file", compare)
    args = ["matrix", command, "--db", str(tmp_path / "project.duckdb"), "--workflow-config", str(settings)]
    if command == "compare":
        args += ["--matrix-file", str(matrix), "--output-file", str(tmp_path / "compare.duckdb")]
    else:
        args += ["--matrix-dir", str(matrices), "--compare-dir", str(tmp_path / "compares")]
    if override:
        args += ["--compare-kernel", override]
    result = CliRunner().invoke(cli.cli, args)
    assert result.exit_code == 0, result.output
    assert calls[0]["compare_kernel"] == expected
    assert calls[0]["calculate"] == "genome_ani"


@pytest.mark.parametrize("override,expected", [(None, "block"), ("pairwise", "anchor")])
def test_slurm_child_gets_resolved_kernel(tmp_path, monkeypatch, override, expected):
    matrices = tmp_path / "matrices"; matrices.mkdir()
    for genome in ("genome_a", "genome_b"):
        (matrices / f"{genome}.h5").write_text("fixture")
    settings = tmp_path / "workflow.toml"
    settings.write_text(
        '[matrix_compare]\ncalculate="genome_ani"\ncompare_kernel="gemm"\n'
        '[stages.matrix_compare]\nexecution="slurm"\nworkers=2\n'
    )
    scripts = []
    def runner(command, **kwargs):
        assert command[0] == "sbatch"
        scripts.append(Path(command[-1]).read_text())
        return subprocess.CompletedProcess(command, 0, stdout="123\n", stderr="")
    monkeypatch.setattr(cli, "WorkflowRuntime", lambda settings, **kwargs: WorkflowRuntime(settings, runner=runner, **kwargs))
    args = ["matrix", "sync-compare", "--db", str(tmp_path / "project.duckdb"),
            "--matrix-dir", str(matrices), "--compare-dir", str(tmp_path / "compares"),
            "--workflow-config", str(settings)]
    if override:
        args += ["--compare-kernel", override]
    result = CliRunner().invoke(cli.cli, args)
    assert result.exit_code == 0, result.output
    assert len(scripts) == 2
    for script in scripts:
        child = shlex.split(script.splitlines()[-1])
        assert child[child.index("--compare-kernel") + 1] == expected
        assert child[child.index("--calculate") + 1] == "genome_ani"
        assert "--no-register" in child


@pytest.mark.parametrize("calculate,ani_method", [
    ("all", "popani"), ("gene_ani", "popani"), ("max_identical_stretch", "popani"),
    ("genome_ani", "cosani_0.5"),
])
def test_invalid_block_combinations_rejected_before_dispatch(tmp_path, monkeypatch, calculate, ani_method):
    settings = tmp_path / "workflow.toml"
    settings.write_text('[stages.matrix_compare]\nexecution="slurm"\n')
    def unexpected(*args, **kwargs):
        raise AssertionError("Invalid kernel combinations must not submit jobs")
    monkeypatch.setattr(cli, "_dispatch_matrix_sync_compare_jobs", unexpected)
    result = CliRunner().invoke(cli.cli, [
        "matrix", "sync-compare", "--db", str(tmp_path / "project.duckdb"),
        "--matrix-dir", str(tmp_path / "matrices"), "--compare-dir", str(tmp_path / "compares"),
        "--workflow-config", str(settings), "--compare-kernel", "block",
        "--calculate", calculate, "--ani-method", ani_method,
    ])
    assert result.exit_code != 0
    assert "block kernel" in result.output
    assert not (tmp_path / "compares").exists()


@pytest.mark.parametrize("route", ["file", "registry"])
def test_python_api_explicit_kernel_overrides_config(tmp_path, monkeypatch, route):
    matrix = tmp_path / "matrix.h5"
    matrix.write_text("fixture")
    monkeypatch.setattr(workflows, "read_matrix_file_context", lambda _: SimpleNamespace(matrix_file=matrix))
    monkeypatch.setattr(db, "get_matrix_store", lambda *args: SimpleNamespace(matrix_file=matrix))
    monkeypatch.setattr(db, "register_matrix_compare", lambda *args, **kwargs: None)
    calls = []
    def compare(**kwargs):
        calls.append(kwargs)
    monkeypatch.setattr(mp, "matrix_compare", compare)
    kwargs = dict(output_file=tmp_path / "result.duckdb", calculate="all", compare_kernel="pairwise",
                  compare_config=config.MatrixCompareConfig(compare_kernel="gemm", target_queue_size=3))
    if route == "file":
        workflows.compare_matrix_file(None, matrix_file=matrix, register=False, **kwargs)
    else:
        workflows.compare_matrix_from_registry(None, matrix_id="genome", **kwargs)
    assert calls[0]["compare_kernel"] == "anchor"
    assert calls[0]["target_queue_size"] == 3


@pytest.mark.parametrize("ani_method,storage_mode", [("popani", "bitmask"), ("conani", "counts")])
def test_real_zipstrain_block_matches_anchor_and_resumes_with_aliases(tmp_path, ani_method, storage_mode):
    profiles = tmp_path / "profiles"; profiles.mkdir()
    for sample, a, c in (("S1", [8, 0, 8, 0], [0, 8, 0, 0]),
                         ("S2", [8, 0, 8, 8], [0, 8, 8, 0]),
                         ("S3", [0, 8, 8, 8], [8, 0, 0, 0])):
        pl.DataFrame({"chrom": ["chr"] * 4, "genome": ["genome"] * 4, "gene": ["gene"] * 4, "pos": [1, 2, 3, 4],
                      "A": a, "C": c, "G": [0] * 4, "T": [0] * 4}).write_parquet(profiles / f"{sample}.parquet")
    bed, stb = tmp_path / "reference.bed", tmp_path / "reference.stb"
    bed.write_text("chr\t0\t4\n"); stb.write_text("chr\tgenome\n")
    matrix = tmp_path / "matrix.h5"
    mp.build_matrix_hdf5(profile_dir=profiles, output_file=matrix, bed_file=bed, stb_file=stb,
                         storage_mode=storage_mode, memory_limit_gb=1)
    results = []
    for kernel, alias in (("anchor", "pairwise"), ("block", "gemm")):
        output = tmp_path / f"{kernel}.duckdb"
        # Old configs and canonical names must address the same completed pairs.
        workflows.compare_matrix_file(
            None, matrix_file=matrix, output_file=output, calculate="ani", ani_method=ani_method,
            compare_kernel=alias, register=False, memory_limit_gb=1, logger=WorkflowLogger(),
        )
        workflows.compare_matrix_file(
            None, matrix_file=matrix, output_file=output, calculate="genome_ani", ani_method=ani_method,
            compare_kernel=kernel, register=False, memory_limit_gb=1, logger=WorkflowLogger(),
        )
        with duckdb.connect(str(output), read_only=True) as conn:
            results.append(conn.execute(
                "SELECT sample_1,sample_2,total_positions,share_allele_pos,genome_ani "
                "FROM matrix_compare_results ORDER BY sample_1,sample_2"
            ).fetchall())
    assert len(results[0]) == 3
    assert results[0] == results[1]

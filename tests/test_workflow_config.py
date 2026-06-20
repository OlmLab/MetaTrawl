from pathlib import Path
import subprocess

import pytest

from metatrawl.config import load_workflow_config
from metatrawl.execution import WorkflowRuntime
from metatrawl.logging import WorkflowLogger


def test_toml_config_controls_stages_and_matrix_compare(tmp_path: Path) -> None:
    path = tmp_path / "workflow.toml"
    path.write_text('''
sample_workers = 9
[stages.bowtie_build]
workers = 1
threads = 12
execution = "slurm"
environment = { OMP_NUM_THREADS = "12" }
[stages.bowtie_build.slurm]
time = "03:00:00"
memory_gb = 48
partition = "compute"
[matrix_compare]
calculate = "ani+ibs"
memory_limit_gb = 48.5
target_queue_size = 3
loader_executor_kind = "process"
''')
    config = load_workflow_config(path, threads=4, sample_count=20)
    assert config.sample_workers == 9
    assert config.stage("bowtie_build").workers == 1
    assert config.stage("bowtie_build").threads == 12
    assert config.stage("bowtie_build").execution == "slurm"
    assert config.stage("bowtie_build").slurm.memory_gb == 48
    assert config.matrix_compare.calculate == "ani+ibs"
    assert config.matrix_compare.memory_limit_gb == 48.5
    assert config.matrix_compare.kwargs() == {"target_queue_size": 3, "loader_executor_kind": "process"}


def test_json_config_is_supported(tmp_path: Path) -> None:
    path = tmp_path / "workflow.json"
    path.write_text('{"sample_workers": 2, "stages": {"alignment": {"workers": 1, "threads": 8}}}')
    config = load_workflow_config(path, threads=4, sample_count=4)
    assert config.sample_workers == 2
    assert config.stage("alignment").workers == 1
    assert config.stage("alignment").threads == 8


def test_config_rejects_unknown_stage(tmp_path: Path) -> None:
    path = tmp_path / "workflow.toml"
    path.write_text("[stages.typo]\nworkers = 2\n")
    with pytest.raises(ValueError, match="Unknown workflow stages: typo"):
        load_workflow_config(path, threads=2, sample_count=2)


def test_config_rejects_nonpositive_matrix_memory_limit(tmp_path: Path) -> None:
    path = tmp_path / "workflow.toml"
    path.write_text("[matrix_compare]\nmemory_limit_gb = 0\n")
    with pytest.raises(ValueError, match="matrix_compare.memory_limit_gb must be a positive number"):
        load_workflow_config(path, threads=2, sample_count=2)


def test_slurm_runtime_writes_script_and_waits(tmp_path: Path) -> None:
    path = tmp_path / "workflow.toml"
    path.write_text('''
[stages.alignment]
execution = "slurm"
workers = 1
threads = 16
environment = { CUDA_VISIBLE_DEVICES = "0" }
[stages.alignment.slurm]
time = "02:00:00"
memory_gb = 64
account = "lab"
''')
    config = load_workflow_config(path, threads=1, sample_count=1)
    calls = []

    def runner(command, **kwargs):
        calls.append((command, kwargs))
        return subprocess.CompletedProcess(command, 0, stdout="123", stderr="")

    runtime = WorkflowRuntime(config, state_dir=tmp_path, logger=WorkflowLogger(), runner=runner)
    runtime.run_shell("alignment", "bowtie2 --version | cat", sample="SRR1")
    command = calls[0][0]
    assert command[:3] == ["sbatch", "--wait", "--parsable"]
    assert "--cpus-per-task" in command and command[command.index("--cpus-per-task") + 1] == "16"
    assert "--account" in command and command[command.index("--account") + 1] == "lab"
    script = Path(command[-1]).read_text()
    assert "export CUDA_VISIBLE_DEVICES=0" in script
    assert "bowtie2 --version | cat" in script

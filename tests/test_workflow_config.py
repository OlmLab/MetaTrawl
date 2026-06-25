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
retries = 2
retry_delay_seconds = 0.5
environment = { OMP_NUM_THREADS = "12" }
[stages.bowtie_build.slurm]
time = "03:00:00"
memory_gb = 48
partition = "compute"
[matrix_compare]
calculate = "ani+ibs"
genome = "GCF_1"
backend = "mps"
memory_limit_gb = 48.5
target_queue_size = 3
loader_executor_kind = "process"
[profile]
min_mapq = 20
min_baseq = 25
min_read_ani = 0.97
read_inclusion = "proper-pairs"
''')
    config = load_workflow_config(path, threads=4, sample_count=20)
    assert config.sample_workers == 9
    assert config.stage("bowtie_build").workers == 1
    assert config.stage("bowtie_build").threads == 12
    assert config.stage("bowtie_build").execution == "slurm"
    assert config.stage("bowtie_build").retries == 2
    assert config.stage("bowtie_build").retry_delay_seconds == 0.5
    assert config.stage("bowtie_build").slurm.memory_gb == 48
    assert config.matrix_compare.calculate == "ani+ibs"
    assert config.matrix_compare.genome == "GCF_1"
    assert config.matrix_compare.backend == "mps"
    assert config.matrix_compare.memory_limit_gb == 48.5
    assert config.matrix_compare.kwargs() == {"target_queue_size": 3, "loader_executor_kind": "process"}
    assert config.profile.min_mapq == 20
    assert config.profile.min_baseq == 25
    assert config.profile.min_read_ani == 0.97
    assert config.profile.read_inclusion == "proper-pairs"


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


def test_config_rejects_negative_stage_retries(tmp_path: Path) -> None:
    path = tmp_path / "workflow.toml"
    path.write_text("[stages.alignment]\nretries = -1\n")
    with pytest.raises(ValueError, match="stages.alignment.retries must be a non-negative integer"):
        load_workflow_config(path, threads=2, sample_count=2)


def test_config_rejects_invalid_profile_read_ani(tmp_path: Path) -> None:
    path = tmp_path / "workflow.toml"
    path.write_text("[profile]\nmin_read_ani = 1.5\n")
    with pytest.raises(ValueError, match="profile.min_read_ani must be between 0 and 1"):
        load_workflow_config(path, threads=2, sample_count=2)


def test_config_rejects_invalid_profile_read_inclusion(tmp_path: Path) -> None:
    path = tmp_path / "workflow.toml"
    path.write_text('[profile]\nread_inclusion = "mapped"\n')
    with pytest.raises(ValueError, match="profile.read_inclusion must be one of"):
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


def test_slurm_runtime_retries_failed_stage(tmp_path: Path) -> None:
    path = tmp_path / "workflow.toml"
    path.write_text('''
[stages.alignment]
execution = "slurm"
workers = 1
threads = 8
retries = 2
retry_delay_seconds = 0
[stages.alignment.slurm]
time = "01:00:00"
memory_gb = 16
''')
    config = load_workflow_config(path, threads=1, sample_count=1)
    calls = []

    def runner(command, **kwargs):
        calls.append((command, kwargs))
        if len(calls) < 3:
            raise subprocess.CalledProcessError(1, command, stderr="PREEMPTED")
        return subprocess.CompletedProcess(command, 0, stdout="123", stderr="")

    runtime = WorkflowRuntime(config, state_dir=tmp_path, logger=WorkflowLogger(), runner=runner)
    runtime.run_shell("alignment", "bowtie2 --version", sample="SRR1")

    assert len(calls) == 3
    assert all(call[0][0] == "sbatch" for call in calls)

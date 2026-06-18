"""Stage-aware local and Slurm command execution."""
from __future__ import annotations
import os
from pathlib import Path
import shlex
import subprocess
import threading
import uuid
from typing import Callable, Sequence

from metatrawl.config import WorkflowConfig
from metatrawl.logging import WorkflowLogger

RunCallable = Callable[..., subprocess.CompletedProcess]

class WorkflowRuntime:
    """Apply per-stage concurrency, resources, and environments to commands."""

    def __init__(self, config: WorkflowConfig, *, state_dir: Path, logger: WorkflowLogger, runner: RunCallable = subprocess.run) -> None:
        self.config = config
        self.state_dir = Path(state_dir)
        self.logger = logger
        self.runner = runner
        self._limits = {name: threading.BoundedSemaphore(stage.workers) for name, stage in config.stages.items()}

    def threads(self, stage: str) -> int:
        return self.config.stage(stage).threads

    def run(self, stage: str, cmd: Sequence[str], *, sample: str, stdout_file: Path | None = None):
        setting = self.config.stage(stage)
        with self._limits[stage]:
            self.logger.emit(sample=sample, step=stage.replace("_", "-"), status="executing", execution=setting.execution, threads=setting.threads)
            if setting.execution == "slurm":
                return self._run_slurm(stage, list(cmd), sample=sample, stdout_file=stdout_file)
            return self._run_local(stage, list(cmd), sample=sample, stdout_file=stdout_file)

    def run_shell(self, stage: str, command: str, *, sample: str) -> None:
        setting = self.config.stage(stage)
        with self._limits[stage]:
            self.logger.emit(sample=sample, step=stage.replace("_", "-"), status="executing", execution=setting.execution, threads=setting.threads)
            if setting.execution == "slurm":
                self._run_slurm(stage, command, sample=sample)
            else:
                self._invoke(command, sample=sample, stage=stage, shell=True, env=self._environment(stage))

    def _run_local(self, stage: str, cmd: list[str], *, sample: str, stdout_file: Path | None) -> None:
        if stdout_file is None:
            return self._invoke(cmd, sample=sample, stage=stage, env=self._environment(stage))
        stdout_file.parent.mkdir(parents=True, exist_ok=True)
        with stdout_file.open("w") as handle:
            return self._invoke(cmd, sample=sample, stage=stage, stdout=handle, env=self._environment(stage))

    def _run_slurm(self, stage: str, command: list[str] | str, *, sample: str, stdout_file: Path | None = None) -> None:
        setting = self.config.stage(stage)
        job_dir = self.state_dir / ".metatrawl_slurm" / stage
        job_dir.mkdir(parents=True, exist_ok=True)
        token = uuid.uuid4().hex[:10]
        script = job_dir / f"{_safe(sample)}-{token}.sh"
        stdout_log = stdout_file or job_dir / f"{_safe(sample)}-{token}.out"
        stderr_log = job_dir / f"{_safe(sample)}-{token}.err"
        command_text = command if isinstance(command, str) else shlex.join(command)
        exports = "\n".join(f"export {key}={shlex.quote(value)}" for key, value in setting.environment.items())
        script.write_text(f"#!/bin/bash\nset -euo pipefail\n{exports}\n{command_text}\n")
        script.chmod(0o700)
        sbatch = [
            "sbatch", "--wait", "--parsable", "--job-name", f"mt-{stage}-{_safe(sample)[:24]}",
            "--cpus-per-task", str(setting.threads), "--mem", f"{setting.slurm.memory_gb}G",
            "--time", setting.slurm.time, "--output", str(stdout_log), "--error", str(stderr_log),
        ]
        if setting.slurm.partition:
            sbatch.extend(["--partition", setting.slurm.partition])
        if setting.slurm.account:
            sbatch.extend(["--account", setting.slurm.account])
        for key, value in setting.slurm.extra.items():
            sbatch.extend([f"--{key.replace('_', '-')}", value])
        sbatch.append(str(script))
        return self._invoke(sbatch, sample=sample, stage=stage, env=os.environ.copy())

    def _environment(self, stage: str) -> dict[str, str]:
        environment = os.environ.copy()
        environment.update(self.config.stage(stage).environment)
        return environment

    def _invoke(self, command, *, sample: str, stage: str, **kwargs) -> None:
        try:
            return self.runner(command, check=True, capture_output="stdout" not in kwargs, text=True, **kwargs)
        except subprocess.CalledProcessError as exc:
            stderr = getattr(exc, "stderr", None) or ""
            rendered = command if isinstance(command, str) else shlex.join(command)
            raise RuntimeError(f"sample={sample} step={stage.replace('_', '-')} command failed: {rendered}\n{stderr}") from exc

def _safe(value: str) -> str:
    return "".join(char if char.isalnum() or char in "._-" else "_" for char in value) or "job"

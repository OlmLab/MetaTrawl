"""Stage-aware local and Slurm command execution."""
from __future__ import annotations
import math
import os
from pathlib import Path
import shlex
import subprocess
import threading
import time
import uuid
from typing import Callable, Sequence

from metatrawl.config import WorkflowConfig
from metatrawl.logging import WorkflowLogger

RunCallable = Callable[..., subprocess.CompletedProcess]


class SlurmJobError(RuntimeError):
    """A failed Slurm job with its best available terminal-state classification."""

    def __init__(self, message: str, *, state: str) -> None:
        super().__init__(message)
        self.state = state

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
            return self._run_with_retries(
                stage,
                sample=sample,
                action=lambda memory_gb, time_limit: (
                    self._run_slurm(
                        stage,
                        list(cmd),
                        sample=sample,
                        stdout_file=stdout_file,
                        memory_gb=memory_gb,
                        time_limit=time_limit,
                    )
                    if setting.execution == "slurm"
                    else self._run_local(stage, list(cmd), sample=sample, stdout_file=stdout_file)
                ),
            )

    def run_shell(self, stage: str, command: str, *, sample: str) -> None:
        setting = self.config.stage(stage)
        with self._limits[stage]:
            self._run_with_retries(
                stage,
                sample=sample,
                action=lambda memory_gb, time_limit: (
                    self._run_slurm(
                        stage,
                        command,
                        sample=sample,
                        memory_gb=memory_gb,
                        time_limit=time_limit,
                    )
                    if setting.execution == "slurm"
                    else self._invoke(command, sample=sample, stage=stage, shell=True, env=self._environment(stage))
                ),
            )

    def _run_with_retries(self, stage: str, *, sample: str, action):
        setting = self.config.stage(stage)
        attempts = setting.retries + 1
        last_error: Exception | None = None
        memory_gb = setting.slurm.memory_gb
        time_limit = setting.slurm.time
        for attempt in range(1, attempts + 1):
            self.logger.emit(
                sample=sample,
                step=stage.replace("_", "-"),
                status="executing",
                execution=setting.execution,
                threads=setting.threads,
                attempt=attempt,
                max_attempts=attempts,
                memory_gb=memory_gb if setting.execution == "slurm" else None,
                time_limit=time_limit if setting.execution == "slurm" else None,
            )
            try:
                return action(memory_gb, time_limit)
            except RuntimeError as exc:
                last_error = exc
                if attempt >= attempts:
                    break
                slurm_state = exc.state if isinstance(exc, SlurmJobError) else None
                if slurm_state == "OUT_OF_MEMORY":
                    memory_gb = max(
                        memory_gb,
                        math.ceil(memory_gb * setting.slurm.memory_retry_coefficient),
                    )
                elif slurm_state == "TIMEOUT":
                    time_limit = _scale_slurm_time(
                        time_limit,
                        setting.slurm.time_retry_coefficient,
                    )
                self.logger.emit(
                    sample=sample,
                    step=stage.replace("_", "-"),
                    status="retrying",
                    execution=setting.execution,
                    attempt=attempt,
                    next_attempt=attempt + 1,
                    max_attempts=attempts,
                    slurm_state=slurm_state,
                    next_memory_gb=memory_gb if setting.execution == "slurm" else None,
                    next_time_limit=time_limit if setting.execution == "slurm" else None,
                    error=exc,
                )
                if setting.retry_delay_seconds > 0:
                    time.sleep(setting.retry_delay_seconds)
        raise last_error if last_error is not None else RuntimeError(f"sample={sample} step={stage} failed")

    def _run_local(self, stage: str, cmd: list[str], *, sample: str, stdout_file: Path | None) -> None:
        if stdout_file is None:
            return self._invoke(cmd, sample=sample, stage=stage, env=self._environment(stage))
        stdout_file.parent.mkdir(parents=True, exist_ok=True)
        with stdout_file.open("w") as handle:
            return self._invoke(cmd, sample=sample, stage=stage, stdout=handle, env=self._environment(stage))

    def _run_slurm(
        self,
        stage: str,
        command: list[str] | str,
        *,
        sample: str,
        stdout_file: Path | None = None,
        memory_gb: int | None = None,
        time_limit: str | None = None,
    ) -> None:
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
            "--cpus-per-task", str(setting.threads), "--mem", f"{memory_gb or setting.slurm.memory_gb}G",
            "--time", time_limit or setting.slurm.time, "--output", str(stdout_log), "--error", str(stderr_log),
        ]
        if setting.slurm.partition:
            sbatch.extend(["--partition", setting.slurm.partition])
        if setting.slurm.account:
            sbatch.extend(["--account", setting.slurm.account])
        for key, value in setting.slurm.extra.items():
            sbatch.extend([f"--{key.replace('_', '-')}", value])
        sbatch.append(str(script))
        try:
            return self.runner(
                sbatch,
                check=True,
                capture_output=True,
                text=True,
                env=os.environ.copy(),
            )
        except subprocess.CalledProcessError as exc:
            scheduler_stderr = getattr(exc, "stderr", None) or ""
            child_stderr = _read_log(stderr_log)
            failure_text = "\n".join(part for part in (scheduler_stderr, child_stderr) if part).strip()
            rendered = shlex.join(sbatch)
            raise SlurmJobError(
                f"sample={sample} step={stage.replace('_', '-')} command failed: {rendered}"
                + (f"\n{failure_text}" if failure_text else ""),
                state=_classify_slurm_failure(failure_text),
            ) from exc

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


def _read_log(path: Path) -> str:
    try:
        return path.read_text(errors="replace")
    except OSError:
        return ""


def _classify_slurm_failure(text: str) -> str:
    normalized = text.upper().replace("-", "_")
    if any(
        marker in normalized
        for marker in ("OUT_OF_MEMORY", "OUT OF MEMORY", "OOM_KILL", "OOM KILL")
    ):
        return "OUT_OF_MEMORY"
    if "TIMEOUT" in normalized or "TIME LIMIT" in normalized:
        return "TIMEOUT"
    if "PREEMPT" in normalized:
        return "PREEMPTED"
    return "FAILED"


def _scale_slurm_time(value: str, coefficient: float) -> str:
    if coefficient == 1:
        return value
    seconds = _parse_slurm_time(value)
    return _format_slurm_time(max(seconds, math.ceil(seconds * coefficient)))


def _parse_slurm_time(value: str) -> int:
    value = str(value).strip()
    if not value:
        raise ValueError("Slurm time cannot be empty.")
    days = 0
    clock = value
    if "-" in value:
        day_text, clock = value.split("-", 1)
        days = int(day_text)
    parts = [int(part) for part in clock.split(":")]
    if len(parts) == 1:
        hours, minutes, seconds = 0, parts[0], 0
    elif len(parts) == 2:
        hours, minutes, seconds = 0, parts[0], parts[1]
    elif len(parts) == 3:
        hours, minutes, seconds = parts
    else:
        raise ValueError(f"Unsupported Slurm time value: {value}")
    if min(days, hours, minutes, seconds) < 0 or minutes >= 60 or seconds >= 60:
        raise ValueError(f"Unsupported Slurm time value: {value}")
    return (((days * 24) + hours) * 60 + minutes) * 60 + seconds


def _format_slurm_time(total_seconds: int) -> str:
    days, remainder = divmod(total_seconds, 24 * 60 * 60)
    hours, remainder = divmod(remainder, 60 * 60)
    minutes, seconds = divmod(remainder, 60)
    clock = f"{hours:02d}:{minutes:02d}:{seconds:02d}"
    return f"{days}-{clock}" if days else clock

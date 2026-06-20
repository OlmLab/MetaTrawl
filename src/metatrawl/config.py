"""Validated workflow configuration for MetaTrawl execution stages."""
from __future__ import annotations
from dataclasses import dataclass, field
import json
from pathlib import Path
import tomllib
from typing import Any

STAGE_NAMES = ("sra_download", "sylph", "genome_download", "prodigal", "prepare_profile", "bowtie_build", "alignment", "profile")

@dataclass(frozen=True)
class SlurmConfig:
    time: str = "01:00:00"
    memory_gb: int = 8
    partition: str | None = None
    account: str | None = None
    extra: dict[str, str] = field(default_factory=dict)

@dataclass(frozen=True)
class StageConfig:
    workers: int = 1
    threads: int = 1
    execution: str = "local"
    environment: dict[str, str] = field(default_factory=dict)
    slurm: SlurmConfig = field(default_factory=SlurmConfig)

@dataclass(frozen=True)
class MatrixCompareConfig:
    calculate: str | None = None
    memory_limit_gb: float | None = None
    anchor_queue_size: int | None = None
    target_queue_size: int | None = None
    result_transfer_batch_size: int | None = None
    loader_executor_kind: str | None = None
    writer_executor_kind: str | None = None

    def kwargs(self) -> dict[str, Any]:
        values = {
            "anchor_queue_size": self.anchor_queue_size,
            "target_queue_size": self.target_queue_size,
            "result_transfer_batch_size": self.result_transfer_batch_size,
            "loader_executor_kind": self.loader_executor_kind,
            "writer_executor_kind": self.writer_executor_kind,
        }
        return {key: value for key, value in values.items() if value is not None}

@dataclass(frozen=True)
class WorkflowConfig:
    sample_workers: int
    stages: dict[str, StageConfig]
    matrix_compare: MatrixCompareConfig = field(default_factory=MatrixCompareConfig)

    def stage(self, name: str) -> StageConfig:
        try:
            return self.stages[name]
        except KeyError as exc:
            raise ValueError(f"Unknown workflow stage: {name}") from exc

    @classmethod
    def legacy(cls, *, threads: int, sample_count: int) -> "WorkflowConfig":
        sample_workers = max(1, min(sample_count or 1, threads))
        stage_threads = max(1, threads // sample_workers)
        stages = {name: StageConfig(workers=sample_workers, threads=stage_threads) for name in STAGE_NAMES}
        stages["genome_download"] = StageConfig(workers=4, threads=1)
        stages["prodigal"] = StageConfig(workers=1, threads=1)
        return cls(sample_workers=sample_workers, stages=stages)

def load_workflow_config(path: Path | None, *, threads: int, sample_count: int) -> WorkflowConfig:
    if path is None:
        return WorkflowConfig.legacy(threads=threads, sample_count=sample_count)
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"Workflow configuration does not exist: {path}")
    with path.open("rb") as handle:
        if path.suffix.lower() == ".toml":
            raw = tomllib.load(handle)
        elif path.suffix.lower() == ".json":
            raw = json.load(handle)
        else:
            raise ValueError("Workflow configuration must use .toml or .json.")
    if not isinstance(raw, dict):
        raise ValueError("Workflow configuration must contain an object/table at its root.")
    return _parse(raw, threads=threads, sample_count=sample_count)

def _parse(raw: dict[str, Any], *, threads: int, sample_count: int) -> WorkflowConfig:
    defaults = WorkflowConfig.legacy(threads=threads, sample_count=sample_count)
    sample_workers = _positive(raw.get("sample_workers", defaults.sample_workers), "sample_workers")
    stage_values = raw.get("stages", {})
    if not isinstance(stage_values, dict):
        raise ValueError("stages must be a TOML table or JSON object.")
    unknown = set(stage_values) - set(STAGE_NAMES)
    if unknown:
        raise ValueError(f"Unknown workflow stages: {', '.join(sorted(unknown))}")
    stages = {}
    for name in STAGE_NAMES:
        default = defaults.stage(name)
        value = stage_values.get(name, {})
        if not isinstance(value, dict):
            raise ValueError(f"stages.{name} must be a table/object.")
        execution = str(value.get("execution", default.execution)).lower()
        if execution not in {"local", "slurm"}:
            raise ValueError(f"stages.{name}.execution must be 'local' or 'slurm'.")
        slurm_value = value.get("slurm", {})
        if not isinstance(slurm_value, dict):
            raise ValueError(f"stages.{name}.slurm must be a table/object.")
        stages[name] = StageConfig(
            workers=_positive(value.get("workers", default.workers), f"stages.{name}.workers"),
            threads=_positive(value.get("threads", default.threads), f"stages.{name}.threads"),
            execution=execution,
            environment=_mapping(value.get("environment", {}), f"stages.{name}.environment"),
            slurm=SlurmConfig(
                time=str(slurm_value.get("time", "01:00:00")),
                memory_gb=_positive(slurm_value.get("memory_gb", 8), f"stages.{name}.slurm.memory_gb"),
                partition=_optional_string(slurm_value.get("partition")),
                account=_optional_string(slurm_value.get("account")),
                extra=_mapping(slurm_value.get("extra", {}), f"stages.{name}.slurm.extra"),
            ),
        )
    compare = raw.get("matrix_compare", {})
    if not isinstance(compare, dict):
        raise ValueError("matrix_compare must be a table/object.")
    return WorkflowConfig(sample_workers=sample_workers, stages=stages, matrix_compare=MatrixCompareConfig(
        calculate=_optional_string(compare.get("calculate")),
        memory_limit_gb=_optional_positive_float(compare.get("memory_limit_gb"), "matrix_compare.memory_limit_gb"),
        anchor_queue_size=_optional_positive(compare.get("anchor_queue_size"), "matrix_compare.anchor_queue_size"),
        target_queue_size=_optional_positive(compare.get("target_queue_size"), "matrix_compare.target_queue_size"),
        result_transfer_batch_size=_optional_positive(compare.get("result_transfer_batch_size"), "matrix_compare.result_transfer_batch_size"),
        loader_executor_kind=_optional_string(compare.get("loader_executor_kind")),
        writer_executor_kind=_optional_string(compare.get("writer_executor_kind")),
    ))

def _positive(value: Any, key: str) -> int:
    try:
        result = int(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{key} must be a positive integer.") from exc
    if result < 1:
        raise ValueError(f"{key} must be a positive integer.")
    return result

def _optional_positive(value: Any, key: str) -> int | None:
    return None if value is None else _positive(value, key)

def _optional_positive_float(value: Any, key: str) -> float | None:
    if value is None:
        return None
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{key} must be a positive number.") from exc
    if result <= 0:
        raise ValueError(f"{key} must be a positive number.")
    return result

def _optional_string(value: Any) -> str | None:
    return None if value in (None, "") else str(value)

def _mapping(value: Any, key: str) -> dict[str, str]:
    if not isinstance(value, dict):
        raise ValueError(f"{key} must be a table/object.")
    return {str(name): str(item) for name, item in value.items()}

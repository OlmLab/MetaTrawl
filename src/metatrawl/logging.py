"""Friendly logging helpers for MetaTrawl workflows."""

from __future__ import annotations

from dataclasses import dataclass, field
import sys
import time


@dataclass
class WorkflowLogger:
    """Emit compact structured logs that stay readable in cluster output."""

    prefix: str = "METATRAWL"
    started_at: float = field(default_factory=time.monotonic)

    def emit(self, *, step: str, status: str, sample: str | None = None, accession: str | None = None, **fields: object) -> None:
        parts = [self.prefix]
        if sample is not None:
            parts.append(f"sample={sample}")
        if accession is not None:
            parts.append(f"accession={accession}")
        parts.append(f"step={step}")
        parts.append(f"status={status}")
        for key, value in fields.items():
            if value is not None:
                parts.append(f"{key}={value}")
        parts.append(f"elapsed={time.monotonic() - self.started_at:.1f}s")
        print(" ".join(parts), file=sys.stderr, flush=True)

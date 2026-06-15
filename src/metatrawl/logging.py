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


class ThrottledMatrixLogger:
    """Emit ZipStrain-style matrix progress without flooding cluster logs."""

    def __init__(self, prefix: str, *, stored_rows: bool = True) -> None:
        self.prefix = prefix
        self.stored_rows = stored_rows
        self.started_at = time.monotonic()
        self.last_percent_bucket = -1
        self.last_log_time = 0.0
        self.logged_advance = False
        self.last_processing_detail: str | None = None

    def __call__(self, event: dict[str, object]) -> None:
        phase = str(event.get("phase", ""))
        completed = int(event.get("completed", 0))
        total = int(event.get("total", 0))
        now = time.monotonic()

        if phase in {"start", "done"}:
            self.last_percent_bucket = -1 if total <= 0 else int((completed / max(total, 1)) * 100) // 5
            self.last_log_time = now
            if phase == "start":
                self.logged_advance = False
                self.last_processing_detail = None
            self._emit(f"{phase.upper()}", event, completed=completed, total=total)
            return

        if phase == "processing":
            detail = event.get("detail")
            if detail is None:
                return
            detail_str = str(detail)
            if detail_str == self.last_processing_detail and (now - self.last_log_time) < 5.0:
                return
            self.last_processing_detail = detail_str
            self.last_log_time = now
            self._emit("PROCESSING", event, completed=completed, total=total, detail=detail_str)
            return

        if phase != "advance":
            return

        percent_bucket = int((completed / max(total, 1)) * 100) // 5 if total > 0 else 0
        should_log = (
            not self.logged_advance
            or total <= 20
            or completed == total
            or percent_bucket > self.last_percent_bucket
            or (now - self.last_log_time) >= 5.0
        )
        if not should_log:
            return

        self.logged_advance = True
        self.last_percent_bucket = percent_bucket
        self.last_log_time = now
        self._emit(
            "PROGRESS",
            event,
            completed=completed,
            total=total,
            percent=f"{(completed / max(total, 1)) * 100:.1f}" if total > 0 else "0.0",
        )

    def _emit(self, label: str, event: dict[str, object], **fields: object) -> None:
        parts = [f"{self.prefix} {label}"]
        parts.append(f"elapsed={time.monotonic() - self.started_at:.1f}s")
        for key, value in fields.items():
            if value is not None:
                parts.append(f"{key}={value}")
        if self.stored_rows:
            parts.append(f"stored_rows={event.get('stored_rows', 0)}")
        genome = event.get("genome")
        if genome:
            parts.append(f"genome={genome}")
        print(" ".join(parts), file=sys.stderr, flush=True)

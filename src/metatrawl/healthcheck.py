"""Dependency health checks for MetaTrawl users."""

from __future__ import annotations

from dataclasses import dataclass
import importlib
from importlib import metadata
import shutil
from typing import Callable

from rich.console import Console
from rich.table import Table


@dataclass(frozen=True)
class Check:
    """One dependency check result."""

    name: str
    status: str
    detail: str


ExecutableProbe = Callable[[str], tuple[bool, str]]
PackageProbe = Callable[[str, str | None], tuple[bool, str]]


def probe_executable(command: str) -> tuple[bool, str]:
    """Return whether an executable is present in PATH."""
    path = shutil.which(command)
    if path is None:
        return False, "not found in PATH"
    return True, path


def probe_package(import_name: str, package_name: str | None = None) -> tuple[bool, str]:
    """Return whether a Python package can be imported."""
    try:
        module = importlib.import_module(import_name)
    except ImportError:
        return False, "not installed"

    try:
        version = metadata.version(package_name or import_name)
    except metadata.PackageNotFoundError:
        version = getattr(module, "__version__", None)
    return True, f"installed ({version})" if version else "installed"


def collect_checks(
    *,
    executable_probe: ExecutableProbe = probe_executable,
    package_probe: PackageProbe = probe_package,
) -> list[Check]:
    """Collect the small dependency report shown by ``metatrawl test``."""
    checks: list[Check] = []
    for command in ("zipstrain", "sylph", "samtools", "bowtie2", "prefetch", "fasterq-dump"):
        ok, detail = executable_probe(command)
        checks.append(Check(name=command, status="ok" if ok else "missing", detail=detail))
    for import_name, package_name in (("torch", "torch"), ("h5py", "h5py")):
        ok, detail = package_probe(import_name, package_name)
        checks.append(Check(name=import_name, status="ok" if ok else "missing", detail=detail))
    return checks


def render_checks(checks: list[Check], *, console: Console | None = None) -> None:
    """Render dependency checks as a compact table."""
    console = console or Console()
    table = Table(title="MetaTrawl health check")
    table.add_column("Dependency", style="cyan")
    table.add_column("Status", style="bold")
    table.add_column("Detail")
    for check in checks:
        status = "[green]OK[/green]" if check.status == "ok" else "[red]MISSING[/red]"
        table.add_row(check.name, status, check.detail)
    console.print(table)
    missing = [check.name for check in checks if check.status != "ok"]
    if missing:
        console.print(f"Missing optional or required tools: {', '.join(missing)}")
    else:
        console.print("All checked dependencies are available.")

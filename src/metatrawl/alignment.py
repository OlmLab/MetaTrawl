"""Alignment adapters with a shared SAM-to-sorted-BAM handoff."""

from __future__ import annotations

from abc import ABC, abstractmethod
from importlib.metadata import entry_points
from pathlib import Path
import shlex


def _nonempty(path: Path) -> bool:
    return path.is_file() and path.stat().st_size > 0


def _require_commands(commands: tuple[str, ...], *, engine: str) -> str:
    """Check the actual execution environment, including a Slurm worker's PATH."""
    checks = []
    for command in commands:
        message = (
            f"Alignment engine '{engine}' requires '{command}' in the executing stage's PATH. "
            "Install it separately or configure the stage environment."
        )
        checks.append(
            f"command -v {shlex.quote(command)} >/dev/null 2>&1 || "
            f"{{ printf '%s\\n' {shlex.quote(message)} >&2; exit 127; }}"
        )
    return "; ".join(checks)


class AlignmentEngine(ABC):
    """Protocol for reference indexing, SAM alignment, and scientific provenance.

    Adapters supply argv, index sidecars, and a scientific contract. MetaTrawl
    owns execution/retries, mapped-read filtering, sorting, BAM validation, and
    publication. Optional adapters can register the ``metatrawl.alignment_engines``
    entry point; no aligner is installed or imported by selecting a built-in.
    """

    name: str

    @abstractmethod
    def index_command(self, reference: Path, threads: int) -> list[str]:
        """Return argv that builds an index alongside the reference."""

    @abstractmethod
    def index_paths(self, reference: Path) -> list[Path]:
        """Return sidecars required for a complete index."""

    @abstractmethod
    def alignment_command(self, reference: Path, reads: list[Path], threads: int) -> list[str]:
        """Return argv that writes alignments as SAM to stdout."""

    @abstractmethod
    def scientific_contract(self) -> dict:
        """Identify alignment semantics, excluding operational thread counts."""

    def index_complete(self, reference: Path) -> bool:
        return all(_nonempty(path) for path in self.index_paths(reference))

    def cleanup_paths(self, reference: Path) -> list[Path]:
        return self.index_paths(reference)

    def prepare_index(self, runtime, reference: Path, *, sample: str) -> None:
        """Build with the existing index-stage resource and retry policy."""
        command = self.index_command(reference, runtime.threads("bowtie_build"))
        script = _require_commands((command[0],), engine=self.name) + "; exec " + shlex.join(command)
        runtime.run_shell("bowtie_build", script, sample=sample)

    def bam_command(self, reference: Path, reads: list[Path], threads: int, bam_file: Path) -> str:
        """Keep only mapped records, sort, and propagate any pipeline failure."""
        if threads < 1 or not reads:
            raise ValueError("Alignment requires positive threads and at least one FASTQ file.")
        align = self.alignment_command(reference, reads, threads)
        script = (
            _require_commands((align[0], "samtools"), engine=self.name)
            + "; " + shlex.join(align)
            + " | samtools view -bS -F 4 -"
            + f" | samtools sort -@ {threads} -o {shlex.quote(str(bam_file))} -"
        )
        return "bash -o pipefail -c " + shlex.quote(script)


class Bowtie2AlignmentEngine(AlignmentEngine):
    """Bowtie2's unchanged sensitive end-to-end defaults."""

    name = "bowtie2"

    def index_command(self, reference: Path, threads: int) -> list[str]:
        return ["bowtie2-build", "--threads", str(threads), str(reference), str(reference)]

    def index_paths(self, reference: Path) -> list[Path]:
        return [
            Path(f"{reference}.{suffix}.{extension}")
            for extension in ("bt2", "bt2l")
            for suffix in ("1", "2", "3", "4", "rev.1", "rev.2")
        ]

    def index_complete(self, reference: Path) -> bool:
        paths = self.index_paths(reference)
        return all(_nonempty(path) for path in paths[:6]) or all(_nonempty(path) for path in paths[6:])

    def prepare_index(self, runtime, reference: Path, *, sample: str) -> None:
        # Preserve the existing Bowtie2 command/runner interface.
        try:
            runtime.run("bowtie_build", self.index_command(reference, runtime.threads("bowtie_build")), sample=sample)
        except FileNotFoundError as exc:
            raise RuntimeError(
                f"sample={sample} step=bowtie-build requires bowtie2-build; install Bowtie2 in the stage environment."
            ) from exc

    def alignment_command(self, reference: Path, reads: list[Path], threads: int) -> list[str]:
        if len(reads) == 2:
            inputs = ["-1", str(reads[0]), "-2", str(reads[1])]
        else:
            inputs = ["-U", ",".join(str(path) for path in reads)]
        return ["bowtie2", "-x", str(reference), *inputs, "--threads", str(threads)]

    def scientific_contract(self) -> dict:
        return {"aligner": self.name, "preset": "default-sensitive-end-to-end", "samtools_exclude_flags": 4}


class BwaMem3AlignmentEngine(AlignmentEngine):
    """Optional user-installed bwa-mem3, using normal index/mem defaults."""

    name = "bwa-mem3"

    def index_command(self, reference: Path, threads: int) -> list[str]:
        return [self.name, "index", "-t", str(threads), str(reference)]

    def index_paths(self, reference: Path) -> list[Path]:
        return [Path(f"{reference}.{suffix}") for suffix in ("amb", "ann", "pac", "bwt.2bit.64")]

    def cleanup_paths(self, reference: Path) -> list[Path]:
        return [*self.index_paths(reference), Path(f"{reference}.0123")]

    def alignment_command(self, reference: Path, reads: list[Path], threads: int) -> list[str]:
        if len(reads) not in (1, 2):
            raise ValueError("bwa-mem3 requires one single-end FASTQ or two paired-end FASTQs per sample.")
        return [self.name, "mem", "-t", str(threads), str(reference), *(str(path) for path in reads)]

    def scientific_contract(self) -> dict:
        return {"aligner": self.name, "preset": "default-mem", "samtools_exclude_flags": 4}


def get_alignment_engine(name: str = "bowtie2") -> AlignmentEngine:
    """Resolve a built-in adapter or an installed adapter entry point."""
    if not isinstance(name, str) or not name.strip() or any(char.isspace() for char in name.strip()):
        raise ValueError("profile.alignment_engine must be a non-empty engine name, for example bowtie2 or bwa-mem3.")
    name = name.strip()
    builtins = {"bowtie2": Bowtie2AlignmentEngine, "bwa-mem3": BwaMem3AlignmentEngine}
    if name in builtins:
        return builtins[name]()
    matches = list(entry_points(group="metatrawl.alignment_engines", name=name))
    if not matches:
        raise ValueError(f"Unknown alignment engine {name!r}. Choose bowtie2 or bwa-mem3, or install its adapter plugin.")
    if len(matches) != 1:
        raise ValueError(f"Multiple alignment adapters registered as {name!r}.")
    engine = matches[0].load()()
    if not isinstance(engine, AlignmentEngine) or engine.name != name:
        raise ValueError(f"Alignment adapter {name!r} must implement AlignmentEngine with the matching name.")
    return engine

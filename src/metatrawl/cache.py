"""Genome cache preparation for concurrent MetaTrawl workers."""

from __future__ import annotations

from concurrent.futures import Future, ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path
import shutil
import subprocess
import threading
from typing import Callable

import polars as pl

from metatrawl.logging import WorkflowLogger


Downloader = Callable[[str, Path], None]
ProdigalRunner = Callable[[Path, Path], None]


@dataclass(frozen=True)
class PreparedReference:
    """Per-sample temporary reference files returned by the cache manager."""

    reference_fasta: Path
    gene_fasta: Path
    stb_file: Path | None = None


class GenomeCache:
    """Single-process authoritative writer for genome/prodigal cache files."""

    def __init__(
        self,
        cache_dir: Path,
        *,
        downloader: Downloader | None = None,
        prodigal_runner: ProdigalRunner | None = None,
        logger: WorkflowLogger | None = None,
    ) -> None:
        self.cache_dir = Path(cache_dir)
        self.genomes_dir = self.cache_dir / "genomes"
        self.genes_dir = self.cache_dir / "genes"
        self.genomes_dir.mkdir(parents=True, exist_ok=True)
        self.genes_dir.mkdir(parents=True, exist_ok=True)
        self.downloader = downloader or download_genome_with_datasets
        self.prodigal_runner = prodigal_runner or run_prodigal_gene_fasta
        self.logger = logger or WorkflowLogger()
        self._executor = ThreadPoolExecutor(max_workers=4)
        self._lock = threading.Lock()
        self._in_flight: dict[str, Future[tuple[Path, Path]]] = {}

    def prepare_accession(self, accession: str) -> tuple[Path, Path]:
        """Ensure one accession has cached genome and gene FASTA files."""
        accession = accession.strip()
        if not accession:
            raise ValueError("Empty accession requested from genome cache.")
        genome_fasta = self.genome_fasta_path(accession)
        gene_fasta = self.gene_fasta_path(accession)
        if genome_fasta.exists() and gene_fasta.exists():
            self.logger.emit(accession=accession, step="cache", status="cached")
            return genome_fasta, gene_fasta
        with self._lock:
            future = self._in_flight.get(accession)
            if future is None:
                self.logger.emit(accession=accession, step="cache", status="queued")
                future = self._executor.submit(self._prepare_accession_uncached, accession)
                self._in_flight[accession] = future
        try:
            return future.result()
        finally:
            if future.done():
                with self._lock:
                    self._in_flight.pop(accession, None)

    def prepare_reference(self, *, accessions: list[str], output_dir: Path, sample: str | None = None) -> PreparedReference:
        """Build per-sample concatenated genome and gene FASTA files."""
        clean_accessions = _dedupe(accessions)
        if not clean_accessions:
            raise ValueError("Cannot prepare a reference with no accessions.")
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        self.logger.emit(sample=sample, step="cache", status="start", accessions=len(clean_accessions))
        prepared = [self.prepare_accession(accession) for accession in clean_accessions]
        reference_fasta = output_dir / "reference.fna"
        gene_fasta = output_dir / "genes.fna"
        stb_file = output_dir / "reference.stb"
        _concatenate_fastas([pair[0] for pair in prepared], reference_fasta)
        _concatenate_fastas([pair[1] for pair in prepared], gene_fasta)
        _write_stb_file(list(zip(clean_accessions, [pair[0] for pair in prepared])), stb_file)
        self.logger.emit(sample=sample, step="cache", status="done", accessions=len(clean_accessions))
        return PreparedReference(reference_fasta=reference_fasta, gene_fasta=gene_fasta, stb_file=stb_file)

    def genome_fasta_path(self, accession: str) -> Path:
        return self.genomes_dir / f"{_safe_name(accession)}.fna"

    def gene_fasta_path(self, accession: str) -> Path:
        return self.genes_dir / f"{_safe_name(accession)}.genes.fna"

    def _prepare_accession_uncached(self, accession: str) -> tuple[Path, Path]:
        genome_fasta = self.genome_fasta_path(accession)
        gene_fasta = self.gene_fasta_path(accession)
        if not genome_fasta.exists():
            self.logger.emit(accession=accession, step="download", status="start")
            tmp_genome = genome_fasta.with_suffix(".tmp.fna")
            self.downloader(accession, tmp_genome)
            _atomic_publish(tmp_genome, genome_fasta)
            self.logger.emit(accession=accession, step="download", status="done", file=genome_fasta)
        if not gene_fasta.exists():
            self.logger.emit(accession=accession, step="prodigal", status="start")
            tmp_gene = gene_fasta.with_suffix(".tmp.fna")
            self.prodigal_runner(genome_fasta, tmp_gene)
            _atomic_publish(tmp_gene, gene_fasta)
            self.logger.emit(accession=accession, step="prodigal", status="done", file=gene_fasta)
        return genome_fasta, gene_fasta


def prepare_cache_reference(
    *,
    cache_dir: Path,
    accessions: list[str],
    output_dir: Path,
    logger: WorkflowLogger | None = None,
) -> PreparedReference:
    """Prepare a per-sample reference using the default cache manager."""
    return GenomeCache(cache_dir, logger=logger).prepare_reference(accessions=accessions, output_dir=output_dir)


def read_accessions_file(accessions_file: Path) -> list[str]:
    """Read accessions from a one-column CSV/TSV/text file."""
    path = Path(accessions_file)
    if path.suffix.lower() == ".csv":
        df = pl.read_csv(path)
        column = "accession" if "accession" in df.columns else df.columns[0]
        return [str(value) for value in df[column].to_list()]
    return [line.strip().split()[0] for line in path.read_text().splitlines() if line.strip()]


def download_genome_with_datasets(accession: str, output_fasta: Path) -> None:
    """Download a genome FASTA with NCBI datasets CLI when available."""
    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    archive = output_fasta.with_suffix(".zip")
    tmp_dir = output_fasta.parent / f"{output_fasta.stem}.datasets"
    try:
        try:
            subprocess.run(
                ["datasets", "download", "genome", "accession", accession, "--filename", str(archive)],
                check=True,
                capture_output=True,
                text=True,
            )
        except OSError as exc:
            raise RuntimeError(
                "Failed to execute NCBI Datasets CLI command `datasets`. "
                "This usually means the `datasets` binary on PATH is for the wrong OS/CPU architecture, "
                "is corrupted, or is not a real executable. Reinstall the NCBI Datasets CLI binary for "
                "your platform and confirm `datasets --version` works."
            ) from exc
        shutil.unpack_archive(str(archive), str(tmp_dir))
        fasta_files = sorted(tmp_dir.rglob("*.fna"))
        if not fasta_files:
            raise RuntimeError(f"datasets did not produce a FASTA for accession {accession}")
        shutil.copy2(fasta_files[0], output_fasta)
    finally:
        if archive.exists():
            archive.unlink()
        if tmp_dir.exists():
            shutil.rmtree(tmp_dir)


def run_prodigal_gene_fasta(genome_fasta: Path, output_gene_fasta: Path) -> None:
    """Run Prodigal and keep only nucleotide gene FASTA output."""
    output_gene_fasta.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        ["prodigal", "-i", str(genome_fasta), "-d", str(output_gene_fasta), "-p", "meta", "-q"],
        check=True,
        capture_output=True,
        text=True,
    )


def _concatenate_fastas(inputs: list[Path], output: Path) -> None:
    tmp = output.with_suffix(output.suffix + ".tmp")
    with tmp.open("w") as out:
        for fasta in inputs:
            out.write(Path(fasta).read_text())
            if not out.tell():
                continue
            out.write("\n")
    _atomic_publish(tmp, output)


def _write_stb_file(accession_fastas: list[tuple[str, Path]], output: Path) -> None:
    tmp = output.with_suffix(output.suffix + ".tmp")
    with tmp.open("w") as handle:
        for accession, fasta in accession_fastas:
            for scaffold in _fasta_headers(fasta):
                handle.write(f"{scaffold}\t{accession}\n")
    _atomic_publish(tmp, output)


def _fasta_headers(fasta: Path) -> list[str]:
    headers: list[str] = []
    with Path(fasta).open() as handle:
        for line in handle:
            if line.startswith(">"):
                headers.append(line[1:].strip().split()[0])
    return headers


def _atomic_publish(tmp: Path, final: Path) -> None:
    final.parent.mkdir(parents=True, exist_ok=True)
    tmp.replace(final)


def _safe_name(value: str) -> str:
    return "".join(char if char.isalnum() or char in "._-" else "_" for char in value)


def _dedupe(values: list[str]) -> list[str]:
    seen: set[str] = set()
    result: list[str] = []
    for value in values:
        clean = value.strip()
        if clean and clean not in seen:
            seen.add(clean)
            result.append(clean)
    return result

"""Genome cache preparation for concurrent MetaTrawl workers."""

from __future__ import annotations

from concurrent.futures import Future, ThreadPoolExecutor
from dataclasses import dataclass
import gzip
from pathlib import Path
import shutil
import subprocess
import threading
import time
from typing import Callable

import polars as pl

from metatrawl.logging import WorkflowLogger


Downloader = Callable[[str, Path], None]
ProdigalRunner = Callable[[Path, Path], None]

DATASETS_RETRY_DELAYS = (5.0, 20.0, 60.0)
PRODIGAL_RETRY_DELAYS = (2.0, 10.0)
DATASETS_PERMANENT_FAILURE_MARKERS = (
    "not found",
    "invalid accession",
    "invalid assembly",
    "suppressed",
    "withdrawn",
    "no assemblies found",
    "no genome found",
    "does not exist",
)


class GenomeUnavailableError(RuntimeError):
    """NCBI does not provide a genome FASTA for the requested accession."""


@dataclass(frozen=True)
class PreparedReference:
    """Per-sample temporary reference files returned by the cache manager."""

    reference_fasta: Path
    gene_fasta: Path
    stb_file: Path | None = None


@dataclass(frozen=True)
class MatrixReferenceFiles:
    """Reusable reference files required to build ZipStrain matrices."""

    reference_fasta: Path
    gene_fasta: Path
    bed_file: Path
    stb_file: Path
    gene_range_table: Path
    accessions: tuple[str, ...]


@dataclass(frozen=True)
class MatrixRequirementSync:
    """Per-accession matrix helper files refreshed from an existing cache."""

    accessions: tuple[str, ...]
    bed_dir: Path
    stb_dir: Path
    gene_range_dir: Path


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
        self.beds_dir = self.cache_dir / "beds"
        self.stb_dir = self.cache_dir / "stb"
        self.gene_ranges_dir = self.cache_dir / "gene_ranges"
        self.genomes_dir.mkdir(parents=True, exist_ok=True)
        self.genes_dir.mkdir(parents=True, exist_ok=True)
        self.beds_dir.mkdir(parents=True, exist_ok=True)
        self.stb_dir.mkdir(parents=True, exist_ok=True)
        self.gene_ranges_dir.mkdir(parents=True, exist_ok=True)
        self.downloader = downloader or download_genome_with_datasets
        self.prodigal_runner = prodigal_runner or run_prodigal_gene_fasta
        self.logger = logger or WorkflowLogger()
        self._executor = ThreadPoolExecutor(max_workers=4)
        self._lock = threading.Lock()
        self._prodigal_lock = threading.Lock()
        self._in_flight: dict[str, Future[tuple[Path, Path]]] = {}

    def prepare_accession(self, accession: str) -> tuple[Path, Path]:
        """Ensure one accession has cached genome and gene FASTA files."""
        accession = accession.strip()
        if not accession:
            raise ValueError("Empty accession requested from genome cache.")
        genome_fasta = self.genome_fasta_path(accession)
        gene_fasta = self.gene_fasta_path(accession)
        bed_file = self.bed_file_path(accession)
        stb_file = self.stb_file_path(accession)
        gene_range_table = self.gene_range_table_path(accession)
        if genome_fasta.exists() and gene_fasta.exists() and bed_file.exists() and stb_file.exists() and gene_range_table.exists():
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
        prepared: list[tuple[str, Path, Path]] = []
        unavailable: list[str] = []
        for accession in clean_accessions:
            try:
                genome_fasta, gene_fasta = self.prepare_accession(accession)
            except GenomeUnavailableError as exc:
                unavailable.append(accession)
                self.logger.emit(
                    sample=sample,
                    accession=accession,
                    step="cache",
                    status="skipped-unavailable",
                    error=exc,
                )
                continue
            prepared.append((accession, genome_fasta, gene_fasta))
        if not prepared:
            raise GenomeUnavailableError(
                "None of the requested accessions produced a usable genome FASTA: "
                + ", ".join(unavailable)
            )
        reference_fasta = output_dir / "reference.fna"
        gene_fasta = output_dir / "genes.fna"
        stb_file = output_dir / "reference.stb"
        _concatenate_fastas([item[1] for item in prepared], reference_fasta)
        _concatenate_fastas([item[2] for item in prepared], gene_fasta)
        _write_stb_file([(item[0], item[1]) for item in prepared], stb_file)
        self.logger.emit(
            sample=sample,
            step="cache",
            status="done",
            accessions=len(prepared),
            unavailable=len(unavailable),
        )
        return PreparedReference(reference_fasta=reference_fasta, gene_fasta=gene_fasta, stb_file=stb_file)

    def genome_fasta_path(self, accession: str) -> Path:
        return self.genomes_dir / f"{_safe_name(accession)}.fna"

    def gene_fasta_path(self, accession: str) -> Path:
        return self.genes_dir / f"{_safe_name(accession)}.genes.fna"

    def bed_file_path(self, accession: str) -> Path:
        return self.beds_dir / f"{_safe_name(accession)}.bed"

    def stb_file_path(self, accession: str) -> Path:
        return self.stb_dir / f"{_safe_name(accession)}.stb"

    def gene_range_table_path(self, accession: str) -> Path:
        return self.gene_ranges_dir / f"{_safe_name(accession)}.gene_ranges.tsv"

    def _prepare_accession_uncached(self, accession: str) -> tuple[Path, Path]:
        genome_fasta = self.genome_fasta_path(accession)
        gene_fasta = self.gene_fasta_path(accession)
        bed_file = self.bed_file_path(accession)
        stb_file = self.stb_file_path(accession)
        gene_range_table = self.gene_range_table_path(accession)
        if not genome_fasta.exists():
            self.logger.emit(accession=accession, step="download", status="start")
            tmp_genome = genome_fasta.with_suffix(".tmp.fna")
            self.downloader(accession, tmp_genome)
            _atomic_publish(tmp_genome, genome_fasta)
            self.logger.emit(accession=accession, step="download", status="done", file=genome_fasta)
        if not gene_fasta.exists():
            self.logger.emit(accession=accession, step="prodigal", status="start")
            tmp_gene = gene_fasta.with_suffix(".tmp.fna")
            # Keep downloads concurrent, but avoid simultaneous Prodigal model training.
            with self._prodigal_lock:
                self.prodigal_runner(genome_fasta, tmp_gene)
            _atomic_publish(tmp_gene, gene_fasta)
            self.logger.emit(accession=accession, step="prodigal", status="done", file=gene_fasta)
        if not bed_file.exists():
            self.logger.emit(accession=accession, step="bed", status="start")
            _write_bed_file([genome_fasta], bed_file, max_interval=500_000)
            self.logger.emit(accession=accession, step="bed", status="done", file=bed_file)
        if not stb_file.exists():
            self.logger.emit(accession=accession, step="stb", status="start")
            _write_stb_file([(accession, genome_fasta)], stb_file)
            self.logger.emit(accession=accession, step="stb", status="done", file=stb_file)
        if not gene_range_table.exists():
            self.logger.emit(accession=accession, step="gene-range", status="start")
            _write_gene_range_table([gene_fasta], gene_range_table)
            self.logger.emit(accession=accession, step="gene-range", status="done", file=gene_range_table)
        return genome_fasta, gene_fasta


def sync_matrix_requirement_files(
    *,
    cache_dir: Path,
    accessions: list[str] | None = None,
    genome: str | None = None,
    max_bed_interval: int = 500_000,
) -> MatrixRequirementSync:
    """Refresh per-genome BED, STB, and gene-range files from cached FASTAs."""
    if genome is not None and accessions is not None:
        raise ValueError("Use either genome or accessions, not both")
    if max_bed_interval <= 0:
        raise ValueError("max_bed_interval must be greater than zero")

    manager = GenomeCache(cache_dir)
    cached = _cached_genome_accessions(manager.genomes_dir)
    selected = _dedupe([genome]) if genome is not None else (_dedupe(accessions) if accessions is not None else cached)
    if not selected:
        raise ValueError(f"No genome FASTA files found in: {manager.genomes_dir}")

    missing_genomes = [accession for accession in selected if not manager.genome_fasta_path(accession).exists()]
    missing_genes = [accession for accession in selected if not manager.gene_fasta_path(accession).exists()]
    if missing_genomes:
        raise FileNotFoundError("Missing genome FASTA for: " + ", ".join(missing_genomes))
    if missing_genes:
        raise FileNotFoundError("Missing Prodigal gene FASTA for: " + ", ".join(missing_genes))

    for accession in selected:
        _write_bed_file([manager.genome_fasta_path(accession)], manager.bed_file_path(accession), max_interval=max_bed_interval)
        _write_stb_file([(accession, manager.genome_fasta_path(accession))], manager.stb_file_path(accession))
        _write_gene_range_table([manager.gene_fasta_path(accession)], manager.gene_range_table_path(accession))

    return MatrixRequirementSync(
        accessions=tuple(selected),
        bed_dir=manager.beds_dir,
        stb_dir=manager.stb_dir,
        gene_range_dir=manager.gene_ranges_dir,
    )


def prepare_cache_reference(
    *,
    cache_dir: Path,
    accessions: list[str],
    output_dir: Path,
    logger: WorkflowLogger | None = None,
) -> PreparedReference:
    """Prepare a per-sample reference using the default cache manager."""
    return GenomeCache(cache_dir, logger=logger).prepare_reference(accessions=accessions, output_dir=output_dir)


def build_matrix_reference_files(
    *,
    genome_dir: Path,
    gene_dir: Path,
    output_dir: Path,
    accessions: list[str] | None = None,
    genome: str | None = None,
    max_bed_interval: int = 500_000,
) -> MatrixReferenceFiles:
    """Build matrix BED/STB/gene-range inputs from cached FASTA directories."""
    genome_dir = Path(genome_dir)
    gene_dir = Path(gene_dir)
    output_dir = Path(output_dir)
    if not genome_dir.is_dir():
        raise FileNotFoundError(f"Genome directory does not exist: {genome_dir}")
    if not gene_dir.is_dir():
        raise FileNotFoundError(f"Gene directory does not exist: {gene_dir}")
    if max_bed_interval <= 0:
        raise ValueError("max_bed_interval must be greater than zero")
    if genome is not None and accessions is not None:
        raise ValueError("Use either genome or accessions, not both")

    genome_files = {
        path.name.removesuffix(".fna"): path
        for path in sorted(genome_dir.glob("*.fna"))
        if not path.name.endswith(".genes.fna")
    }
    gene_files = {
        path.name.removesuffix(".genes.fna"): path
        for path in sorted(gene_dir.glob("*.genes.fna"))
    }
    if genome is not None:
        selected = _dedupe([genome])
    else:
        selected = _dedupe(accessions) if accessions is not None else sorted(genome_files)
    if not selected:
        raise ValueError(f"No genome FASTA files found in: {genome_dir}")

    missing_genomes = [accession for accession in selected if accession not in genome_files]
    missing_genes = [accession for accession in selected if accession not in gene_files]
    if missing_genomes:
        raise FileNotFoundError("Missing genome FASTA for: " + ", ".join(missing_genomes))
    if missing_genes:
        raise FileNotFoundError("Missing Prodigal gene FASTA for: " + ", ".join(missing_genes))

    output_dir.mkdir(parents=True, exist_ok=True)
    reference_fasta = output_dir / "reference.fna"
    gene_fasta = output_dir / "genes.fna"
    bed_file = output_dir / "genomes_bed_file.bed"
    stb_file = output_dir / "reference.stb"
    gene_range_table = output_dir / "gene_range_table.tsv"

    selected_genomes = [genome_files[accession] for accession in selected]
    selected_genes = [gene_files[accession] for accession in selected]
    _concatenate_fastas(selected_genomes, reference_fasta)
    _concatenate_fastas(selected_genes, gene_fasta)
    _write_stb_file(list(zip(selected, selected_genomes)), stb_file)
    _write_bed_file(selected_genomes, bed_file, max_interval=max_bed_interval)
    _write_gene_range_table(selected_genes, gene_range_table)
    return MatrixReferenceFiles(
        reference_fasta=reference_fasta,
        gene_fasta=gene_fasta,
        bed_file=bed_file,
        stb_file=stb_file,
        gene_range_table=gene_range_table,
        accessions=tuple(selected),
    )


def read_accessions_file(accessions_file: Path) -> list[str]:
    """Read accessions from a one-column CSV/TSV/text file."""
    path = Path(accessions_file)
    if path.suffix.lower() == ".csv":
        df = pl.read_csv(path)
        column = "accession" if "accession" in df.columns else df.columns[0]
        return [str(value) for value in df[column].to_list()]
    return [line.strip().split()[0] for line in path.read_text().splitlines() if line.strip()]


def _cached_genome_accessions(genome_dir: Path) -> list[str]:
    return [
        path.name.removesuffix(".fna")
        for path in sorted(Path(genome_dir).glob("*.fna"))
        if not path.name.endswith(".genes.fna")
    ]


def download_genome_with_datasets(
    accession: str,
    output_fasta: Path,
    *,
    retry_delays: tuple[float, ...] = DATASETS_RETRY_DELAYS,
) -> None:
    """Download a genome FASTA, retrying transient NCBI or network failures."""
    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    archive = output_fasta.with_suffix(".zip")
    tmp_dir = output_fasta.parent / f"{output_fasta.stem}.datasets"
    attempts = len(retry_delays) + 1
    last_error = ""
    for attempt in range(1, attempts + 1):
        _remove_datasets_download_artifacts(archive, tmp_dir)
        if output_fasta.exists():
            output_fasta.unlink()
        try:
            result = subprocess.run(
                [
                    "datasets",
                    "download",
                    "genome",
                    "accession",
                    accession,
                    "--include",
                    "genome",
                    "--filename",
                    str(archive),
                ],
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
        except subprocess.CalledProcessError as exc:
            last_error = _subprocess_detail(exc)
            if _is_permanent_datasets_failure(last_error):
                _remove_datasets_download_artifacts(archive, tmp_dir)
                raise GenomeUnavailableError(
                    f"NCBI Datasets cannot provide accession {accession}: {last_error}"
                ) from exc
            if attempt < attempts:
                time.sleep(retry_delays[attempt - 1])
                continue
            _remove_datasets_download_artifacts(archive, tmp_dir)
            raise RuntimeError(
                f"NCBI Datasets download failed for accession {accession} after {attempts} attempts. "
                f"Last error: {last_error or 'unknown transient failure'}"
            ) from exc

        try:
            if not archive.is_file():
                raise ValueError("datasets command completed without creating an archive")
            shutil.unpack_archive(str(archive), str(tmp_dir))
            fasta_files = sorted(tmp_dir.rglob("*.fna"))
            compressed_fasta_files = sorted(tmp_dir.rglob("*.fna.gz"))
            if fasta_files:
                shutil.copy2(fasta_files[0], output_fasta)
            elif compressed_fasta_files:
                with gzip.open(compressed_fasta_files[0], "rb") as source, output_fasta.open("wb") as destination:
                    shutil.copyfileobj(source, destination)
            else:
                detail = _command_output(result)
                if _is_permanent_datasets_failure(detail):
                    raise GenomeUnavailableError(
                        f"NCBI Datasets has no genome FASTA for accession {accession}: "
                        f"{detail or 'assembly unavailable'}"
                    )
                raise ValueError("downloaded archive contained no genome FASTA")
        except GenomeUnavailableError:
            _remove_datasets_download_artifacts(archive, tmp_dir)
            raise
        except (OSError, shutil.ReadError, ValueError) as exc:
            last_error = str(exc)
            if attempt < attempts:
                time.sleep(retry_delays[attempt - 1])
                continue
            _remove_datasets_download_artifacts(archive, tmp_dir)
            raise GenomeUnavailableError(
                f"NCBI Datasets produced no usable FASTA for accession {accession} "
                f"after {attempts} attempts; skipping this accession. Last error: {last_error}"
            ) from exc
        else:
            _remove_datasets_download_artifacts(archive, tmp_dir)
            return


def _remove_datasets_download_artifacts(archive: Path, tmp_dir: Path) -> None:
    if archive.exists():
        archive.unlink()
    if tmp_dir.exists():
        shutil.rmtree(tmp_dir)


def _command_output(result: subprocess.CompletedProcess[str]) -> str:
    return "\n".join(part.strip() for part in (result.stderr, result.stdout) if part and part.strip())


def _subprocess_detail(exc: subprocess.CalledProcessError) -> str:
    detail = "\n".join(
        part.strip()
        for part in (getattr(exc, "stderr", None), getattr(exc, "stdout", None))
        if part and part.strip()
    )
    return detail or f"datasets exited with status {exc.returncode}"


def _is_permanent_datasets_failure(detail: str) -> bool:
    normalized = detail.lower()
    return any(marker in normalized for marker in DATASETS_PERMANENT_FAILURE_MARKERS)


def run_prodigal_gene_fasta(
    genome_fasta: Path,
    output_gene_fasta: Path,
    *,
    retry_delays: tuple[float, ...] = PRODIGAL_RETRY_DELAYS,
) -> None:
    """Annotate one assembled genome, retrying transient process crashes."""
    output_gene_fasta.parent.mkdir(parents=True, exist_ok=True)
    command = ["prodigal", "-i", str(genome_fasta), "-d", str(output_gene_fasta), "-p", "single", "-q"]
    attempts = len(retry_delays) + 1
    for attempt in range(1, attempts + 1):
        if output_gene_fasta.exists():
            output_gene_fasta.unlink()
        try:
            subprocess.run(command, check=True, capture_output=True, text=True)
            return
        except OSError:
            raise
        except subprocess.CalledProcessError as exc:
            if attempt < attempts:
                time.sleep(retry_delays[attempt - 1])
                continue
            detail = _subprocess_detail(exc)
            signal_detail = f" signal={-exc.returncode}" if exc.returncode < 0 else ""
            raise RuntimeError(
                f"Prodigal failed for {genome_fasta} after {attempts} attempts;"
                f"{signal_detail} last_error={detail}"
            ) from exc


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


def _write_bed_file(inputs: list[Path], output: Path, *, max_interval: int) -> None:
    tmp = output.with_suffix(output.suffix + ".tmp")
    with tmp.open("w") as handle:
        for fasta in inputs:
            for scaffold, length in _fasta_lengths(fasta):
                for start in range(0, length, max_interval):
                    handle.write(f"{scaffold}\t{start}\t{min(start + max_interval, length)}\n")
    _atomic_publish(tmp, output)


def _write_gene_range_table(inputs: list[Path], output: Path) -> None:
    tmp = output.with_suffix(output.suffix + ".tmp")
    with tmp.open("w") as handle:
        for fasta in inputs:
            with fasta.open() as fasta_handle:
                for line in fasta_handle:
                    if not line.startswith(">"):
                        continue
                    parts = line[1:].strip().split()
                    if len(parts) < 5 or parts[1] != "#" or parts[3] != "#":
                        raise ValueError(f"Invalid Prodigal gene FASTA header in {fasta}: {line.strip()}")
                    gene = parts[0]
                    scaffold_parts = gene.rsplit("_", 1)
                    if len(scaffold_parts) != 2:
                        raise ValueError(f"Cannot infer scaffold from Prodigal gene ID: {gene}")
                    handle.write(f"{gene}\t{scaffold_parts[0]}\t{int(parts[2])}\t{int(parts[4])}\n")
    _atomic_publish(tmp, output)


def _fasta_lengths(fasta: Path) -> list[tuple[str, int]]:
    records: list[tuple[str, int]] = []
    scaffold: str | None = None
    length = 0
    with Path(fasta).open() as handle:
        for line in handle:
            if line.startswith(">"):
                if scaffold is not None:
                    records.append((scaffold, length))
                scaffold = line[1:].strip().split()[0]
                length = 0
            else:
                length += len(line.strip())
    if scaffold is not None:
        records.append((scaffold, length))
    return records


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

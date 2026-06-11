from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
import csv
from pathlib import Path
import subprocess
import sys
import time
import types
import zipfile

from click.testing import CliRunner
import duckdb
import polars as pl
import pytest

from metatrawl import cache
from metatrawl import cli
from metatrawl import db as registry
from metatrawl import workflows
from metatrawl.logging import WorkflowLogger


def _matrix_contract_files(tmp_path: Path) -> tuple[Path, Path]:
    bed_file = tmp_path / "genomes.bed"
    stb_file = tmp_path / "genomes.stb"
    bed_file.write_text("contigA\t1\t10\tgenome_a\n")
    stb_file.write_text("contigA\tgenome_a\n")
    return bed_file, stb_file


def _write_bundle_files(tmp_path: Path, run_id: str, *, coverage: float = 2.0, breadth: float = 0.9, ber: float = 0.8, abundance: float = 0.1) -> registry.ProfileBundle:
    profile_file = tmp_path / f"{run_id}.parquet"
    genome_stats_file = tmp_path / f"{run_id}.genome_stats.parquet"
    gene_stats_file = tmp_path / f"{run_id}.gene_stats.parquet"
    sylph_file = tmp_path / f"{run_id}.sylph.csv"

    pl.DataFrame(
        {
            "chrom": ["contigA", "contigA"],
            "genome": ["genome_a", "genome_a"],
            "pos": [1, 2],
            "gene": ["gene1", "gene1"],
            "A": [5.0, 0.0],
            "C": [0.0, 3.0],
            "G": [0.0, 0.0],
            "T": [0.0, 0.0],
        }
    ).write_parquet(profile_file)
    pl.DataFrame(
        {
            "genome": ["genome_a"],
            "coverage": [coverage],
            "breadth": [breadth],
            "ber": [ber],
        }
    ).write_parquet(genome_stats_file)
    pl.DataFrame(
        {
            "genome": ["genome_a"],
            "gene": ["gene1"],
            "coverage": [coverage],
            "breadth": [breadth],
            "ber": [ber],
        }
    ).write_parquet(gene_stats_file)
    pl.DataFrame(
        {
            "genome": ["genome_a"],
            "accession": ["GCF_1"],
            "abundance": [abundance],
        }
    ).write_csv(sylph_file)
    return registry.ProfileBundle(
        run_id=run_id,
        profile_file=profile_file,
        genome_stats_file=genome_stats_file,
        gene_stats_file=gene_stats_file,
        sylph_abundance_file=sylph_file,
    )


def _write_manifest(path: Path, bundles: list[registry.ProfileBundle], *, include_gene: bool = True) -> Path:
    fields = ["run_id", "profile_file", "genome_stats_file", "sylph_abundance_file"]
    if include_gene:
        fields.insert(3, "gene_stats_file")
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for bundle in bundles:
            row = {
                "run_id": bundle.run_id,
                "profile_file": str(bundle.profile_file),
                "genome_stats_file": str(bundle.genome_stats_file),
                "sylph_abundance_file": str(bundle.sylph_abundance_file),
            }
            if include_gene:
                row["gene_stats_file"] = str(bundle.gene_stats_file or "")
            writer.writerow(row)
    return path


def _import_bundle(runner: CliRunner, db_file: Path, bundle: registry.ProfileBundle, *, add_run: bool = False) -> None:
    args = [
        "profiles",
        "import",
        "--db",
        str(db_file),
        "--run-id",
        bundle.run_id,
        "--profile-file",
        str(bundle.profile_file),
        "--genome-stats-file",
        str(bundle.genome_stats_file),
        "--sylph-abundance-file",
        str(bundle.sylph_abundance_file),
    ]
    if bundle.gene_stats_file is not None:
        args.extend(["--gene-stats-file", str(bundle.gene_stats_file)])
    if add_run:
        args.append("--add-run")
    result = runner.invoke(cli.cli, args)
    assert result.exit_code == 0, result.output


def test_init_creates_duckdb_schema(tmp_path: Path) -> None:
    runner = CliRunner()
    db_file = tmp_path / "metatrawl.duckdb"

    result = runner.invoke(cli.cli, ["init", "--db", str(db_file)])

    assert result.exit_code == 0, result.output
    tables = {row[0] for row in duckdb.connect(str(db_file)).execute("SHOW TABLES").fetchall()}
    assert {
        "sra_runs",
        "samples",
        "profiles",
        "profile_positions",
        "genome_stats",
        "gene_stats",
        "sylph_abundance",
        "cache_genomes",
        "matrix_stores",
        "matrix_store_samples",
        "matrix_compares",
    } <= tables


def test_adding_run_ids_is_idempotent(tmp_path: Path) -> None:
    runner = CliRunner()
    db_file = tmp_path / "metatrawl.duckdb"

    first = runner.invoke(cli.cli, ["runs", "add", "--db", str(db_file), "SRR1", "SRR1", "SRR2"])
    second = runner.invoke(cli.cli, ["runs", "add", "--db", str(db_file), "SRR1", "SRR2"])

    assert first.exit_code == 0, first.output
    assert "added=2" in first.output
    assert second.exit_code == 0, second.output
    assert "added=0" in second.output


def test_deleted_run_is_excluded_from_remaining_profiles(tmp_path: Path) -> None:
    runner = CliRunner()
    db_file = tmp_path / "metatrawl.duckdb"

    assert runner.invoke(cli.cli, ["runs", "add", "--db", str(db_file), "SRR1", "SRR2"]).exit_code == 0
    assert runner.invoke(cli.cli, ["runs", "delete", "--db", str(db_file), "SRR2"]).exit_code == 0
    output = tmp_path / "remaining.csv"
    result = runner.invoke(cli.cli, ["profiles", "remaining", "--db", str(db_file), "--output-file", str(output)])

    assert result.exit_code == 0, result.output
    assert output.read_text().splitlines() == ["run_id", "SRR1"]


def test_profiles_import_stores_profile_stats_gene_and_abundance_tables(tmp_path: Path) -> None:
    runner = CliRunner()
    db_file = tmp_path / "metatrawl.duckdb"
    bundle = _write_bundle_files(tmp_path, "SRR1")

    _import_bundle(runner, db_file, bundle, add_run=True)

    with duckdb.connect(str(db_file)) as conn:
        assert conn.execute("SELECT count(*) FROM profile_positions WHERE sample_id = 'SRR1'").fetchone()[0] == 2
        assert conn.execute("SELECT coverage, breadth, ber FROM genome_stats WHERE sample_id = 'SRR1'").fetchone() == (2.0, 0.9, 0.8)
        assert conn.execute("SELECT gene FROM gene_stats WHERE sample_id = 'SRR1'").fetchone() == ("gene1",)
        assert conn.execute("SELECT accession, abundance FROM sylph_abundance WHERE sample_id = 'SRR1'").fetchone() == ("GCF_1", 0.1)
        assert conn.execute("SELECT status FROM samples WHERE sample_id = 'SRR1'").fetchone() == ("complete",)
        count_types = dict(
            conn.execute(
                """
                SELECT column_name, data_type
                FROM information_schema.columns
                WHERE table_name = 'profile_positions'
                  AND column_name IN ('A', 'C', 'G', 'T')
                """
            ).fetchall()
        )
        assert count_types == {
            "A": "USMALLINT",
            "C": "USMALLINT",
            "G": "USMALLINT",
            "T": "USMALLINT",
        }


def test_profiles_import_normalizes_sylph_genome_paths(tmp_path: Path) -> None:
    runner = CliRunner()
    db_file = tmp_path / "metatrawl.duckdb"
    bundle = _write_bundle_files(tmp_path, "SRR1")
    pl.DataFrame(
        {
            "Genome_file": [
                "gtdb_genomes_reps_r232/database/GCF/901/875/305/GCF_901875305.1_genomic.fna.gz"
            ],
            "Taxonomic_abundance": [1.5],
        }
    ).write_csv(bundle.sylph_abundance_file)

    _import_bundle(runner, db_file, bundle, add_run=True)

    with duckdb.connect(str(db_file)) as conn:
        assert conn.execute(
            "SELECT genome, accession FROM sylph_abundance WHERE sample_id = 'SRR1'"
        ).fetchone() == ("GCF_901875305.1", "GCF_901875305.1")


def test_profiles_remaining_after_duckdb_import(tmp_path: Path) -> None:
    runner = CliRunner()
    db_file = tmp_path / "metatrawl.duckdb"
    assert runner.invoke(cli.cli, ["runs", "add", "--db", str(db_file), "SRR1", "SRR2"]).exit_code == 0
    _import_bundle(runner, db_file, _write_bundle_files(tmp_path, "SRR1"))

    output = tmp_path / "remaining.csv"
    result = runner.invoke(cli.cli, ["profiles", "remaining", "--db", str(db_file), "--output-file", str(output)])

    assert result.exit_code == 0, result.output
    assert output.read_text().splitlines() == ["run_id", "SRR2"]


def test_profiles_remaining_reports_all_complete_and_writes_header(tmp_path: Path) -> None:
    runner = CliRunner()
    db_file = tmp_path / "metatrawl.duckdb"
    assert runner.invoke(cli.cli, ["runs", "add", "--db", str(db_file), "SRR1"]).exit_code == 0
    _import_bundle(runner, db_file, _write_bundle_files(tmp_path, "SRR1"))

    output = tmp_path / "remaining.csv"
    result = runner.invoke(cli.cli, ["profiles", "remaining", "--db", str(db_file), "--output-file", str(output)])

    assert result.exit_code == 0, result.output
    assert "All added runs have complete profiles" in result.output
    assert output.read_text().splitlines() == ["run_id"]


def test_profiles_add_rejects_missing_required_columns(tmp_path: Path) -> None:
    runner = CliRunner()
    db_file = tmp_path / "metatrawl.duckdb"
    manifest = tmp_path / "bad_manifest.csv"
    manifest.write_text("run_id,profile_file\nSRR1,profile.parquet\n")

    result = runner.invoke(cli.cli, ["profiles", "add", "--db", str(db_file), "--manifest", str(manifest)])

    assert result.exit_code != 0
    assert "missing required columns" in result.output


def test_profiles_add_rejects_unknown_runs_unless_add_runs_is_used(tmp_path: Path) -> None:
    runner = CliRunner()
    db_file = tmp_path / "metatrawl.duckdb"
    manifest = _write_manifest(tmp_path / "manifest.csv", [_write_bundle_files(tmp_path, "SRR1")])

    rejected = runner.invoke(cli.cli, ["profiles", "add", "--db", str(db_file), "--manifest", str(manifest)])
    accepted = runner.invoke(cli.cli, ["profiles", "add", "--db", str(db_file), "--manifest", str(manifest), "--add-runs"])

    assert rejected.exit_code != 0
    assert "unknown run_id" in rejected.output
    assert accepted.exit_code == 0, accepted.output
    assert "imported=1" in accepted.output


def test_cache_prepare_reuses_cached_genome_and_gene_fasta(tmp_path: Path) -> None:
    runner = CliRunner()
    cache_dir = tmp_path / "cache"
    genomes_dir = cache_dir / "genomes"
    genes_dir = cache_dir / "genes"
    genomes_dir.mkdir(parents=True)
    genes_dir.mkdir(parents=True)
    (genomes_dir / "GCF_1.fna").write_text(">g1\nACGT\n")
    (genes_dir / "GCF_1.genes.fna").write_text(">gene1\nAC\n")
    accessions = tmp_path / "accessions.csv"
    accessions.write_text("accession\nGCF_1\n")

    result = runner.invoke(
        cli.cli,
        ["cache", "prepare", "--cache-dir", str(cache_dir), "--accessions", str(accessions), "--output-dir", str(tmp_path / "ref")],
    )

    assert result.exit_code == 0, result.output
    assert ">g1" in (tmp_path / "ref" / "reference.fna").read_text()
    assert ">gene1" in (tmp_path / "ref" / "genes.fna").read_text()
    assert (tmp_path / "ref" / "reference.stb").read_text() == "g1\tGCF_1\n"


def test_build_matrix_reference_files_from_genome_and_gene_directories(tmp_path: Path) -> None:
    genomes = tmp_path / "genomes"
    genes = tmp_path / "genes"
    genomes.mkdir()
    genes.mkdir()
    (genomes / "GCF_1.fna").write_text(">contig_1\nACGTAC\n")
    (genomes / "GCF_2.fna").write_text(">contig_2\nAAAA\n")
    (genes / "GCF_1.genes.fna").write_text(">contig_1_1 # 2 # 5 # 1 # ID=1\nCGTA\n")
    (genes / "GCF_2.genes.fna").write_text(">contig_2_1 # 1 # 3 # 1 # ID=1\nAAA\n")

    files = cache.build_matrix_reference_files(
        genome_dir=genomes,
        gene_dir=genes,
        output_dir=tmp_path / "matrix_reference",
        accessions=["GCF_2"],
    )

    assert files.accessions == ("GCF_2",)
    assert files.reference_fasta.read_text().startswith(">contig_2")
    assert files.stb_file.read_text() == "contig_2\tGCF_2\n"
    assert files.bed_file.read_text() == "contig_2\t0\t4\n"
    assert files.gene_range_table.read_text() == "contig_2_1\tcontig_2\t1\t3\n"


def test_cache_build_matrix_files_cli_uses_all_cached_accessions(tmp_path: Path) -> None:
    genomes = tmp_path / "genomes"
    genes = tmp_path / "genes"
    genomes.mkdir()
    genes.mkdir()
    (genomes / "GCF_1.fna").write_text(">contig_1\nACGT\n")
    (genes / "GCF_1.genes.fna").write_text(">contig_1_1 # 1 # 4 # 1 # ID=1\nACGT\n")
    output = tmp_path / "matrix_reference"

    result = CliRunner().invoke(
        cli.cli,
        [
            "cache",
            "build-matrix-files",
            "--genome-dir",
            str(genomes),
            "--gene-dir",
            str(genes),
            "--output-dir",
            str(output),
        ],
    )

    assert result.exit_code == 0, result.output
    assert '"accessions": 1' in result.output
    assert (output / "genomes_bed_file.bed").exists()
    assert (output / "reference.stb").exists()
    assert (output / "gene_range_table.tsv").exists()


def test_cache_build_matrix_files_cli_can_select_one_genome(tmp_path: Path) -> None:
    genomes = tmp_path / "genomes"
    genes = tmp_path / "genes"
    genomes.mkdir()
    genes.mkdir()
    for accession in ("GCF_1", "GCF_2"):
        (genomes / f"{accession}.fna").write_text(f">contig_{accession[-1]}\nACGT\n")
        (genes / f"{accession}.genes.fna").write_text(
            f">contig_{accession[-1]}_1 # 1 # 4 # 1 # ID=1\nACGT\n"
        )
    output = tmp_path / "matrix_reference"

    result = CliRunner().invoke(
        cli.cli,
        [
            "cache",
            "build-matrix-files",
            "--genome-dir",
            str(genomes),
            "--gene-dir",
            str(genes),
            "--genome",
            "GCF_2",
            "--output-dir",
            str(output),
        ],
    )

    assert result.exit_code == 0, result.output
    assert '"accessions": 1' in result.output
    assert (output / "reference.stb").read_text() == "contig_2\tGCF_2\n"
    assert ">contig_2" in (output / "reference.fna").read_text()
    assert ">contig_1" not in (output / "reference.fna").read_text()


def test_cache_build_matrix_files_rejects_genome_and_accession_list(tmp_path: Path) -> None:
    accessions = tmp_path / "accessions.txt"
    accessions.write_text("GCF_1\n")

    result = CliRunner().invoke(
        cli.cli,
        [
            "cache",
            "build-matrix-files",
            "--genome-dir",
            str(tmp_path),
            "--gene-dir",
            str(tmp_path),
            "--output-dir",
            str(tmp_path / "output"),
            "--genome",
            "GCF_1",
            "--accessions",
            str(accessions),
        ],
    )

    assert result.exit_code != 0
    assert "either --genome or --accessions" in result.output


def test_datasets_exec_format_error_is_actionable(tmp_path: Path, monkeypatch) -> None:
    def fake_run(*args, **kwargs):
        raise OSError(8, "Exec format error")

    monkeypatch.setattr(cache.subprocess, "run", fake_run)

    with pytest.raises(RuntimeError) as exc_info:
        cache.download_genome_with_datasets("GCF_1", tmp_path / "GCF_1.fna")

    assert "NCBI Datasets CLI" in str(exc_info.value)
    assert "wrong OS/CPU architecture" in str(exc_info.value)
    assert "datasets --version" in str(exc_info.value)


def test_datasets_retries_transient_failure_then_downloads_fasta(tmp_path: Path, monkeypatch) -> None:
    calls = 0

    def fake_run(cmd, **kwargs):
        nonlocal calls
        calls += 1
        if calls < 3:
            raise subprocess.CalledProcessError(
                1,
                cmd,
                stderr="connection reset by peer",
            )
        archive = Path(cmd[cmd.index("--filename") + 1])
        with zipfile.ZipFile(archive, "w") as handle:
            handle.writestr("ncbi_dataset/data/GCF_1/GCF_1_genomic.fna", ">contig\nACGT\n")
        return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

    monkeypatch.setattr(cache.subprocess, "run", fake_run)

    output = tmp_path / "GCF_1.fna"
    cache.download_genome_with_datasets("GCF_1", output, retry_delays=(0, 0))

    assert calls == 3
    assert output.read_text() == ">contig\nACGT\n"


def test_datasets_does_not_retry_permanent_accession_failure(tmp_path: Path, monkeypatch) -> None:
    calls = 0

    def fake_run(cmd, **kwargs):
        nonlocal calls
        calls += 1
        raise subprocess.CalledProcessError(
            1,
            cmd,
            stderr="assembly accession not found",
        )

    monkeypatch.setattr(cache.subprocess, "run", fake_run)

    with pytest.raises(cache.GenomeUnavailableError, match="not found"):
        cache.download_genome_with_datasets(
            "GCA_902373375.1",
            tmp_path / "genome.fna",
            retry_delays=(0, 0, 0),
        )

    assert calls == 1


def test_datasets_retries_ambiguous_empty_archives(tmp_path: Path, monkeypatch) -> None:
    calls = 0

    def fake_run(cmd, **kwargs):
        nonlocal calls
        calls += 1
        archive = Path(cmd[cmd.index("--filename") + 1])
        with zipfile.ZipFile(archive, "w") as handle:
            handle.writestr("README.md", "No genome payload")
        return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

    monkeypatch.setattr(cache.subprocess, "run", fake_run)

    with pytest.raises(RuntimeError, match="after 3 attempts"):
        cache.download_genome_with_datasets(
            "GCA_UNKNOWN",
            tmp_path / "genome.fna",
            retry_delays=(0, 0),
        )

    assert calls == 3


def test_datasets_accepts_compressed_genome_fasta(tmp_path: Path, monkeypatch) -> None:
    import gzip

    def fake_run(cmd, **kwargs):
        archive = Path(cmd[cmd.index("--filename") + 1])
        with zipfile.ZipFile(archive, "w") as handle:
            handle.writestr(
                "ncbi_dataset/data/GCF_1/GCF_1_genomic.fna.gz",
                gzip.compress(b">contig\nACGT\n"),
            )
        return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

    monkeypatch.setattr(cache.subprocess, "run", fake_run)

    output = tmp_path / "GCF_1.fna"
    cache.download_genome_with_datasets("GCF_1", output, retry_delays=())

    assert output.read_text() == ">contig\nACGT\n"


def test_duplicate_accession_requests_share_one_preparation_job(tmp_path: Path) -> None:
    calls = {"download": 0, "prodigal": 0}

    def downloader(accession: str, output: Path) -> None:
        calls["download"] += 1
        time.sleep(0.05)
        output.write_text(f">{accession}\nACGT\n")

    def prodigal(genome: Path, output: Path) -> None:
        calls["prodigal"] += 1
        time.sleep(0.05)
        output.write_text(">gene\nAC\n")

    manager = cache.GenomeCache(tmp_path / "cache", downloader=downloader, prodigal_runner=prodigal)
    with ThreadPoolExecutor(max_workers=2) as executor:
        results = list(executor.map(manager.prepare_accession, ["GCF_1", "GCF_1"]))

    assert results[0] == results[1]
    assert calls == {"download": 1, "prodigal": 1}


def test_prodigal_runner_only_requests_nucleotide_gene_fasta(tmp_path: Path, monkeypatch) -> None:
    calls: list[list[str]] = []

    def fake_run(cmd, check, capture_output, text):
        calls.append(cmd)
        Path(cmd[cmd.index("-d") + 1]).write_text(">gene\nAC\n")
        return subprocess.CompletedProcess(cmd, 0)

    monkeypatch.setattr(cache.subprocess, "run", fake_run)

    cache.run_prodigal_gene_fasta(tmp_path / "genome.fna", tmp_path / "genes.fna")

    assert calls
    assert "-d" in calls[0]
    assert "-a" not in calls[0]
    assert (tmp_path / "genes.fna").exists()
    assert not (tmp_path / "genes.faa").exists()


def test_profile_sra_deletes_temporary_reference_files_after_success(tmp_path: Path, monkeypatch) -> None:
    scratch_dir = tmp_path / "scratch"

    def fake_run(cmd: list[str], *, sample: str, step: str) -> None:
        if cmd[0] == "fasterq-dump":
            sample_dir = scratch_dir / sample
            sample_dir.mkdir(parents=True, exist_ok=True)
            (sample_dir / "accessions.txt").write_text("GCF_1\n")

    class FakeGenomeCache:
        def __init__(self, cache_dir: Path, logger: WorkflowLogger | None = None) -> None:
            self.cache_dir = cache_dir
            self.logger = logger

        def prepare_reference(self, *, accessions: list[str], output_dir: Path, sample: str | None = None) -> cache.PreparedReference:
            output_dir.mkdir(parents=True, exist_ok=True)
            reference = output_dir / "reference.fna"
            genes = output_dir / "genes.fna"
            reference.write_text(">g\nACGT\n")
            genes.write_text(">gene\nAC\n")
            return cache.PreparedReference(reference_fasta=reference, gene_fasta=genes)

    monkeypatch.setattr(workflows, "_run", fake_run)
    monkeypatch.setattr(workflows.cache, "GenomeCache", FakeGenomeCache)

    workflows.profile_sra_runs(
        run_ids=["SRR1"],
        db_file=tmp_path / "metatrawl.duckdb",
        cache_dir=tmp_path / "cache",
        scratch_dir=scratch_dir,
        logger=WorkflowLogger(),
    )

    assert not (scratch_dir / "SRR1").exists()


def test_profile_sra_uses_per_sample_accessions_dir(tmp_path: Path, monkeypatch) -> None:
    scratch_dir = tmp_path / "scratch"
    accessions_dir = tmp_path / "accessions"
    accessions_dir.mkdir()
    (accessions_dir / "SRR1.accessions.txt").write_text("GCF_1\n")
    prepared_accessions: list[list[str]] = []

    def fake_run(cmd: list[str], *, sample: str, step: str) -> None:
        return None

    class FakeGenomeCache:
        def __init__(self, cache_dir: Path, logger: WorkflowLogger | None = None) -> None:
            self.cache_dir = cache_dir
            self.logger = logger

        def prepare_reference(self, *, accessions: list[str], output_dir: Path, sample: str | None = None) -> cache.PreparedReference:
            prepared_accessions.append(accessions)
            output_dir.mkdir(parents=True, exist_ok=True)
            reference = output_dir / "reference.fna"
            genes = output_dir / "genes.fna"
            reference.write_text(">g\nACGT\n")
            genes.write_text(">gene\nAC\n")
            return cache.PreparedReference(reference_fasta=reference, gene_fasta=genes)

    monkeypatch.setattr(workflows, "_run", fake_run)
    monkeypatch.setattr(workflows.cache, "GenomeCache", FakeGenomeCache)

    workflows.profile_sra_runs(
        run_ids=["SRR1"],
        db_file=tmp_path / "metatrawl.duckdb",
        cache_dir=tmp_path / "cache",
        scratch_dir=scratch_dir,
        accessions_dir=accessions_dir,
        logger=WorkflowLogger(),
    )

    assert prepared_accessions == [["GCF_1"]]
    assert not (scratch_dir / "SRR1").exists()


def test_profile_sra_runs_sylph_and_extracts_accessions(tmp_path: Path, monkeypatch) -> None:
    scratch_dir = tmp_path / "scratch"
    output_dir = tmp_path / "outputs"
    sylph_db = tmp_path / "gtdb.syldb"
    sylph_db.write_text("fake syldb")
    prepared_accessions: list[list[str]] = []

    def fake_subprocess_run(cmd, check, capture_output=False, text=False, stdout=None, stderr=None):
        if cmd[0] == "prefetch":
            return subprocess.CompletedProcess(cmd, 0)
        if cmd[0] == "fasterq-dump":
            sample_dir = scratch_dir / "SRR1"
            sample_dir.mkdir(parents=True, exist_ok=True)
            (sample_dir / "SRR1_1.fastq").write_text("@r1\nACGT\n+\n!!!!\n")
            (sample_dir / "SRR1_2.fastq").write_text("@r2\nTGCA\n+\n!!!!\n")
            return subprocess.CompletedProcess(cmd, 0)
        if cmd[:2] == ["sylph", "profile"]:
            stdout.write("Genome_file\tTaxonomic_abundance\n/path/GCF_000001.1_genomic.fna\t1.5\n/path/GCF_000002.1_genomic.fna\t0\n")
            return subprocess.CompletedProcess(cmd, 0)
        raise AssertionError(f"unexpected command: {cmd}")

    class FakeGenomeCache:
        def __init__(self, cache_dir: Path, logger: WorkflowLogger | None = None) -> None:
            self.cache_dir = cache_dir
            self.logger = logger

        def prepare_reference(self, *, accessions: list[str], output_dir: Path, sample: str | None = None) -> cache.PreparedReference:
            prepared_accessions.append(accessions)
            output_dir.mkdir(parents=True, exist_ok=True)
            reference = output_dir / "reference.fna"
            genes = output_dir / "genes.fna"
            reference.write_text(">g\nACGT\n")
            genes.write_text(">gene\nAC\n")
            return cache.PreparedReference(reference_fasta=reference, gene_fasta=genes)

    monkeypatch.setattr(workflows.subprocess, "run", fake_subprocess_run)
    monkeypatch.setattr(workflows.cache, "GenomeCache", FakeGenomeCache)
    monkeypatch.setattr(workflows, "_run_alignment_and_profile", lambda **kwargs: None)

    workflows.profile_sra_runs(
        run_ids=["SRR1"],
        db_file=tmp_path / "metatrawl.duckdb",
        cache_dir=tmp_path / "cache",
        scratch_dir=scratch_dir,
        sylph_db=sylph_db,
        output_dir=output_dir,
        logger=WorkflowLogger(),
    )

    assert prepared_accessions == [["GCF_000001.1"]]
    assert (output_dir / "SRR1.sylph.tsv").read_text().startswith("Genome_file")
    assert not (scratch_dir / "SRR1").exists()


def test_alignment_and_profile_stage_publishes_outputs(tmp_path: Path, monkeypatch) -> None:
    sample_scratch = tmp_path / "scratch" / "SRR1"
    sample_scratch.mkdir(parents=True)
    (sample_scratch / "SRR1_1.fastq").write_text("@r1\nACGT\n+\n!!!!\n")
    (sample_scratch / "SRR1_2.fastq").write_text("@r2\nTGCA\n+\n!!!!\n")
    reference_dir = sample_scratch / "reference"
    reference_dir.mkdir()
    reference = cache.PreparedReference(
        reference_fasta=reference_dir / "reference.fna",
        gene_fasta=reference_dir / "genes.fna",
        stb_file=reference_dir / "reference.stb",
    )
    reference.reference_fasta.write_text(">contig\nACGT\n")
    reference.gene_fasta.write_text(">gene\nAC\n")
    reference.stb_file.write_text("contig\tGCF_1\n")
    run_calls: list[list[str]] = []
    shell_calls: list[str] = []

    def fake_run(cmd: list[str], *, sample: str, step: str) -> None:
        run_calls.append(cmd)
        if cmd[:3] == ["zipstrain", "utilities", "prepare_profiling"]:
            output_dir = Path(cmd[cmd.index("--output-dir") + 1])
            output_dir.mkdir(parents=True, exist_ok=True)
            (output_dir / "genomes_bed_file.bed").write_text("contig\t1\t4\n")
            (output_dir / "gene_range_table.tsv").write_text("gene\tcontig\t1\t2\n")
            (output_dir / "null_model.parquet").write_text("null")
            (output_dir / "profiling_contract.json").write_text("{}")
        if cmd[:3] == ["zipstrain", "utilities", "profile-single"]:
            output_dir = Path(cmd[cmd.index("--output-dir") + 1])
            (output_dir / "SRR1_profile.parquet").write_text("profile")
            (output_dir / "SRR1_genome_stats.parquet").write_text("genome")
            (output_dir / "SRR1_gene_stats.parquet").write_text("gene")

    def fake_shell(*, command: str, sample: str) -> None:
        shell_calls.append(command)

    monkeypatch.setattr(workflows, "_run", fake_run)
    monkeypatch.setattr(workflows, "_run_alignment_shell", fake_shell)

    workflows._run_alignment_and_profile(
        run_id="SRR1",
        sample_scratch=sample_scratch,
        reference=reference,
        output_dir=tmp_path / "outputs",
        threads=2,
        logger=WorkflowLogger(),
    )

    assert any(call[:2] == ["bowtie2-build", "--threads"] for call in run_calls)
    profile_calls = [call for call in run_calls if call[:3] == ["zipstrain", "utilities", "profile-single"]]
    assert profile_calls
    assert "--profiling-contract" not in profile_calls[0]
    assert shell_calls
    assert "bowtie2 -x" in shell_calls[0]
    assert "samtools sort" in shell_calls[0]
    assert (tmp_path / "outputs" / "SRR1.profile.parquet").read_text() == "profile"
    assert (tmp_path / "outputs" / "SRR1.genome_stats.parquet").read_text() == "genome"
    assert (tmp_path / "outputs" / "SRR1.gene_stats.parquet").read_text() == "gene"


def test_alignment_profile_builds_null_model_when_prepare_does_not(tmp_path: Path, monkeypatch) -> None:
    sample_scratch = tmp_path / "scratch" / "SRR1"
    sample_scratch.mkdir(parents=True)
    (sample_scratch / "SRR1.fastq").write_text("@r1\nACGT\n+\n!!!!\n")
    reference_dir = sample_scratch / "reference"
    reference_dir.mkdir()
    reference = cache.PreparedReference(
        reference_fasta=reference_dir / "reference.fna",
        gene_fasta=reference_dir / "genes.fna",
        stb_file=reference_dir / "reference.stb",
    )
    reference.reference_fasta.write_text(">contig\nACGT\n")
    reference.gene_fasta.write_text(">gene\nAC\n")
    reference.stb_file.write_text("contig\tGCF_1\n")
    run_calls: list[list[str]] = []

    def fake_run(cmd: list[str], *, sample: str, step: str) -> None:
        run_calls.append(cmd)
        if cmd[:3] == ["zipstrain", "utilities", "prepare_profiling"]:
            output_dir = Path(cmd[cmd.index("--output-dir") + 1])
            output_dir.mkdir(parents=True, exist_ok=True)
            (output_dir / "genomes_bed_file.bed").write_text("contig\t1\t4\n")
            (output_dir / "gene_range_table.tsv").write_text("gene\tcontig\t1\t2\n")
        if cmd[:3] == ["zipstrain", "utilities", "build-null-model"]:
            Path(cmd[cmd.index("--output-file") + 1]).write_text("null")
        if cmd[:3] == ["zipstrain", "utilities", "profile-single"]:
            output_dir = Path(cmd[cmd.index("--output-dir") + 1])
            (output_dir / "SRR1_profile.parquet").write_text("profile")
            (output_dir / "SRR1_genome_stats.parquet").write_text("genome")
            (output_dir / "SRR1_gene_stats.parquet").write_text("gene")

    monkeypatch.setattr(workflows, "_run", fake_run)
    monkeypatch.setattr(workflows, "_run_alignment_shell", lambda **kwargs: None)

    workflows._run_alignment_and_profile(
        run_id="SRR1",
        sample_scratch=sample_scratch,
        reference=reference,
        output_dir=tmp_path / "outputs",
        threads=1,
        logger=WorkflowLogger(),
    )

    assert any(call[:3] == ["zipstrain", "utilities", "build-null-model"] for call in run_calls)
    profile_call = next(call for call in run_calls if call[:3] == ["zipstrain", "utilities", "profile-single"])
    assert (sample_scratch / "zipstrain_profile" / "null_model.parquet").exists()
    assert profile_call[profile_call.index("--null-model") + 1].endswith("null_model.parquet")


def test_sylph_single_end_reads_are_positional_not_dash_u(tmp_path: Path, monkeypatch) -> None:
    sample_scratch = tmp_path / "scratch" / "SRR1"
    sample_scratch.mkdir(parents=True)
    (sample_scratch / "SRR1.fastq").write_text("@r1\nACGT\n+\n!!!!\n")
    sylph_db = tmp_path / "gtdb.syldb"
    sylph_db.write_text("fake syldb")
    calls: list[list[str]] = []

    def fake_run(cmd, check, stdout, stderr, text):
        calls.append(cmd)
        stdout.write("Genome_file\tTaxonomic_abundance\n/path/GCF_000001.1_genomic.fna\t1.5\n")
        return subprocess.CompletedProcess(cmd, 0)

    monkeypatch.setattr(workflows.subprocess, "run", fake_run)

    workflows._run_sylph_profile(
        sylph_db=sylph_db,
        sample_scratch=sample_scratch,
        run_id="SRR1",
        threads=1,
        output_file=tmp_path / "sylph.tsv",
    )

    assert calls == [["sylph", "profile", str(sylph_db), str(sample_scratch / "SRR1.fastq"), "-t", "1"]]
    assert "-U" not in calls[0]


def test_sync_profiles_remaining_runs_and_imports_outputs(tmp_path: Path, monkeypatch) -> None:
    runner = CliRunner()
    db_file = tmp_path / "metatrawl.duckdb"
    output_dir = tmp_path / "outputs"
    assert runner.invoke(cli.cli, ["runs", "add", "--db", str(db_file), "SRR1", "SRR2"]).exit_code == 0
    _import_bundle(runner, db_file, _write_bundle_files(tmp_path, "SRR1"))

    def fake_profile_sra_runs(**kwargs) -> None:
        assert kwargs["run_ids"] == ["SRR2"]
        output_dir.mkdir(parents=True, exist_ok=True)
        bundle = _write_bundle_files(output_dir, "SRR2")
        bundle.profile_file.rename(output_dir / "SRR2.profile.parquet")

    monkeypatch.setattr(workflows, "profile_sra_runs", fake_profile_sra_runs)

    result = runner.invoke(
        cli.cli,
        [
            "sync",
            "--db",
            str(db_file),
            "--cache-dir",
            str(tmp_path / "cache"),
            "--scratch-dir",
            str(tmp_path / "scratch"),
            "--output-dir",
            str(output_dir),
            "--skip-dependency-check",
        ],
    )

    assert result.exit_code == 0, result.output
    assert "requested=1" in result.output
    assert "imported=1" in result.output
    assert "cleaned_files=4" in result.output
    with duckdb.connect(str(db_file)) as conn:
        assert conn.execute("SELECT sample_id FROM samples ORDER BY sample_id").fetchall() == [("SRR1",), ("SRR2",)]
        assert conn.execute("SELECT count(*) FROM profile_positions WHERE sample_id = 'SRR2'").fetchone()[0] == 2
    assert not (output_dir / "SRR2.profile.parquet").exists()
    assert not (output_dir / "SRR2.genome_stats.parquet").exists()
    assert not (output_dir / "SRR2.gene_stats.parquet").exists()
    assert not (output_dir / "SRR2.sylph.csv").exists()


def test_sync_can_keep_profile_outputs_for_debugging(tmp_path: Path, monkeypatch) -> None:
    runner = CliRunner()
    db_file = tmp_path / "metatrawl.duckdb"
    output_dir = tmp_path / "outputs"
    assert runner.invoke(cli.cli, ["runs", "add", "--db", str(db_file), "SRR1"]).exit_code == 0

    def fake_profile_sra_runs(**kwargs) -> None:
        output_dir.mkdir(parents=True, exist_ok=True)
        bundle = _write_bundle_files(output_dir, "SRR1")
        bundle.profile_file.rename(output_dir / "SRR1.profile.parquet")

    monkeypatch.setattr(workflows, "profile_sra_runs", fake_profile_sra_runs)

    result = runner.invoke(
        cli.cli,
        [
            "sync",
            "--db",
            str(db_file),
            "--cache-dir",
            str(tmp_path / "cache"),
            "--scratch-dir",
            str(tmp_path / "scratch"),
            "--output-dir",
            str(output_dir),
            "--skip-dependency-check",
            "--keep-profile-outputs",
        ],
    )

    assert result.exit_code == 0, result.output
    assert "cleaned_files=0" in result.output
    assert (output_dir / "SRR1.profile.parquet").exists()
    assert (output_dir / "SRR1.genome_stats.parquet").exists()
    assert (output_dir / "SRR1.gene_stats.parquet").exists()
    assert (output_dir / "SRR1.sylph.csv").exists()


def test_matrix_build_filters_samples_exports_selected_and_passes_sparse(tmp_path: Path, monkeypatch) -> None:
    runner = CliRunner()
    db_file = tmp_path / "metatrawl.duckdb"
    assert runner.invoke(cli.cli, ["runs", "add", "--db", str(db_file), "SRR_PASS", "SRR_FAIL"]).exit_code == 0
    _import_bundle(runner, db_file, _write_bundle_files(tmp_path, "SRR_PASS", coverage=5, breadth=0.95, ber=0.9, abundance=0.2))
    _import_bundle(runner, db_file, _write_bundle_files(tmp_path, "SRR_FAIL", coverage=0.1, breadth=0.2, ber=0.1, abundance=0.001))
    bed_file, stb_file = _matrix_contract_files(tmp_path)
    calls: list[dict[str, object]] = []

    def build_matrix_hdf5(**kwargs):
        profile_dir = Path(kwargs["profile_dir"])
        calls.append({**kwargs, "staged": sorted(path.name for path in profile_dir.glob("*.parquet"))})
        Path(kwargs["output_file"]).write_text("matrix")

    zipstrain_module = types.ModuleType("zipstrain")
    matrix_pairs_module = types.ModuleType("zipstrain.matrix_pairs")
    matrix_pairs_module.build_matrix_hdf5 = build_matrix_hdf5
    monkeypatch.setitem(sys.modules, "zipstrain", zipstrain_module)
    monkeypatch.setitem(sys.modules, "zipstrain.matrix_pairs", matrix_pairs_module)
    monkeypatch.setattr(zipstrain_module, "matrix_pairs", matrix_pairs_module, raising=False)

    result = runner.invoke(
        cli.cli,
        [
            "matrix",
            "build",
            "--db",
            str(db_file),
            "--genome",
            "genome_a",
            "--bed-file",
            str(bed_file),
            "--stb-file",
            str(stb_file),
            "--output-file",
            str(tmp_path / "matrix.h5"),
            "--min-coverage",
            "1",
            "--min-breadth",
            "0.5",
            "--min-ber",
            "0.5",
            "--min-sylph-abundance",
            "0.01",
            "--sparse",
        ],
    )

    assert result.exit_code == 0, result.output
    assert "profiles=1" in result.output
    assert calls[0]["staged"] == ["SRR_PASS.parquet"]
    assert calls[0]["sparse"] is True
    with duckdb.connect(str(db_file)) as conn:
        assert conn.execute("SELECT storage_layout, profile_count FROM matrix_stores").fetchall() == [("sparse", 1)]
        assert conn.execute("SELECT sample_id FROM matrix_store_samples").fetchall() == [("SRR_PASS",)]


def test_matrix_build_fails_when_no_complete_profiles_exist(tmp_path: Path) -> None:
    runner = CliRunner()
    db_file = tmp_path / "metatrawl.duckdb"
    assert runner.invoke(cli.cli, ["runs", "add", "--db", str(db_file), "SRR1"]).exit_code == 0
    bed_file, stb_file = _matrix_contract_files(tmp_path)

    result = runner.invoke(
        cli.cli,
        [
            "matrix",
            "build",
            "--db",
            str(db_file),
            "--genome",
            "genome_a",
            "--bed-file",
            str(bed_file),
            "--stb-file",
            str(stb_file),
            "--output-file",
            str(tmp_path / "matrix.h5"),
        ],
    )

    assert result.exit_code != 0
    assert "No complete samples" in result.output


def test_matrix_compare_uses_registered_matrix_store(tmp_path: Path, monkeypatch) -> None:
    runner = CliRunner()
    db_file = tmp_path / "metatrawl.duckdb"
    matrix_file = tmp_path / "matrix.h5"
    matrix_file.write_text("matrix")
    calls: list[dict[str, object]] = []

    with registry.connect(db_file) as conn:
        registry.register_matrix_store(
            conn,
            matrix_id="genome_a",
            genome="genome_a",
            matrix_file=matrix_file,
            profile_count=2,
        )

    def matrix_compare(**kwargs):
        calls.append(kwargs)
        Path(kwargs["output_file"]).write_text("compare")

    zipstrain_module = types.ModuleType("zipstrain")
    matrix_pairs_module = types.ModuleType("zipstrain.matrix_pairs")
    matrix_pairs_module.matrix_compare = matrix_compare
    monkeypatch.setitem(sys.modules, "zipstrain", zipstrain_module)
    monkeypatch.setitem(sys.modules, "zipstrain.matrix_pairs", matrix_pairs_module)
    monkeypatch.setattr(zipstrain_module, "matrix_pairs", matrix_pairs_module, raising=False)

    result = runner.invoke(
        cli.cli,
        [
            "matrix",
            "compare",
            "--db",
            str(db_file),
            "--matrix-id",
            "genome_a",
            "--output-file",
            str(tmp_path / "compare.duckdb"),
            "--calculate",
            "all",
        ],
    )

    assert result.exit_code == 0, result.output
    assert "compare_id=compare" in result.output
    assert calls[0]["matrix_db_file"] == matrix_file
    with duckdb.connect(str(db_file)) as conn:
        assert conn.execute("SELECT compare_id, matrix_id, calculate FROM matrix_compares").fetchall() == [("compare", "genome_a", "all")]


def test_workflow_logs_include_clear_step_context(capsys) -> None:
    logger = WorkflowLogger()

    logger.emit(sample="SRR1", accession="GCF_1", step="cache", status="done", genomes=2)

    captured = capsys.readouterr()
    assert "METATRAWL" in captured.err
    assert "sample=SRR1" in captured.err
    assert "accession=GCF_1" in captured.err
    assert "step=cache" in captured.err
    assert "status=done" in captured.err

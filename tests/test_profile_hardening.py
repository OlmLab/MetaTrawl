"""Regression fixtures for real-world download, parsing, and bundle failures."""

from dataclasses import asdict, replace
import json
import os
from pathlib import Path
import subprocess
import sys

from click.testing import CliRunner
import polars as pl
import pytest

from metatrawl import cli, db, provenance, workflows
from metatrawl.config import NullModelConfig, ProfileConfig, load_workflow_config
from metatrawl.logging import WorkflowLogger
from metatrawl.sylph import read_sylph_table
from test_provenance import bundle


@pytest.mark.parametrize("extension", ["tsv", "csv", "parquet"])
def test_mixed_sylph_fields_and_repeated_headers_roundtrip(tmp_path, extension, capsys):
    path = tmp_path / f"sample.{extension}"
    table = pl.DataFrame({
        "Genome_file": ["GCF_1.1.fna", "Genome_file", "GCA_2.1.fna", "GCF_3.1.fna"],
        "Taxonomic_abundance": ["2.5", "Taxonomic_abundance", "0.0", "7.5"],
        "Eff_lambda": ["0.2", "Eff_lambda", "LOW", "HIGH"],
    })
    if extension == "parquet":
        table.write_parquet(path)
    else:
        table.write_csv(path, separator="\t" if extension == "tsv" else ",")
    before = path.read_bytes()
    assert read_sylph_table(path).height == 3
    assert workflows.extract_accessions_from_sylph_table(path) == ["GCF_1.1", "GCF_3.1"]
    with db.connect(tmp_path / "project.duckdb") as conn:
        db._insert_sylph_abundance(conn, sample_id="SRR1", abundance_file=path)
        assert conn.execute("SELECT accession,abundance FROM sylph_abundance ORDER BY accession").fetchall() == [
            ("GCA_2.1", 0.0), ("GCF_1.1", 2.5), ("GCF_3.1", 7.5),
        ]
    assert path.read_bytes() == before
    assert "repeated-header-skipped" in capsys.readouterr().err


def test_high_after_schema_inference_window_is_valid(tmp_path):
    path = tmp_path / "sylph.tsv"
    path.write_text("Genome_file\tTaxonomic_abundance\tEff_lambda\n" +
                    "GCF_1.1\t1\t0.1\n" * 105 + "GCA_2.1\t1\tHIGH\n")
    assert workflows.extract_accessions_from_sylph_table(path) == ["GCA_2.1", "GCF_1.1"]


@pytest.mark.parametrize("value", ["BROKEN", "Taxonomic_abundance", "NaN", "inf", "-1"])
def test_malformed_abundance_is_not_ignored_as_a_header(tmp_path, value):
    path = tmp_path / "sylph.tsv"
    path.write_text(f"Genome_file\tTaxonomic_abundance\tEff_lambda\nGCF_1.1\t{value}\tHIGH\n")
    with pytest.raises(ValueError, match="Sylph column"):
        workflows.extract_accessions_from_sylph_table(path)
    with db.connect(tmp_path / "project.duckdb") as conn:
        with pytest.raises(ValueError, match="Sylph column"):
            db._insert_sylph_abundance(conn, sample_id="SRR1", abundance_file=path)
        assert conn.execute("SELECT count(*) FROM sylph_abundance").fetchone()[0] == 0


def test_header_only_sylph_has_no_accessions(tmp_path):
    path = tmp_path / "sylph.tsv"
    path.write_text("Genome_file\tTaxonomic_abundance\tEff_lambda\n")
    assert workflows.extract_accessions_from_sylph_table(path) == []


@pytest.mark.parametrize("key,value", [
    ("error_rate", 0), ("error_rate", 1), ("error_rate", float("nan")),
    ("error_rate", True), ("p_threshold", float("inf")), ("p_threshold", -0.1),
    ("p_threshold", 1), ("max_total_reads", 0), ("max_total_reads", 1.5),
    ("max_total_reads", True), ("model_type", "binomial"),
])
def test_null_model_rejects_invalid_parameters(key, value):
    with pytest.raises(ValueError, match="profile.null_model"):
        NullModelConfig(**{key: value})


@pytest.mark.parametrize("extension", ["toml", "json"])
def test_null_model_config_and_provenance(tmp_path, extension):
    path = tmp_path / f"workflow.{extension}"
    values = dict(error_rate=0.02, p_threshold=0.001, max_total_reads=500, model_type="poisson")
    path.write_text(
        '[profile.null_model]\nerror_rate=0.02\np_threshold=0.001\nmax_total_reads=500\nmodel_type="poisson"\n'
        if extension == "toml" else json.dumps({"profile": {"null_model": values}})
    )
    config = load_workflow_config(path, threads=2, sample_count=1)
    assert asdict(config.profile.null_model) == values
    contract = provenance.profiling_contract(config.profile)
    assert contract["null_model"] == values
    assert "null_model" not in contract["profile"]
    assert provenance.null_model_options(config.profile) == [
        "--error-rate", "0.02", "--max-total-reads", "500",
        "--p-threshold", "0.001", "--model-type", "poisson",
    ]


def test_null_model_defaults_preserve_old_contract_fingerprint():
    config = ProfileConfig()
    historical_profile = asdict(config)
    for key in ("backend", "alignment_engine", "null_model"):
        historical_profile.pop(key)
    contract = provenance.profiling_contract(config)
    assert contract["profile"] == historical_profile
    assert contract["null_model"] == provenance.zipstrain_defaults()["null_model"]
    explicit = replace(config, null_model=NullModelConfig(**contract["null_model"]))
    assert provenance.fingerprint(provenance.profiling_contract(explicit)) == provenance.fingerprint(contract)


def test_changed_null_model_rejected_before_import(tmp_path):
    with db.connect(tmp_path / "project.duckdb") as conn:
        db.add_runs(conn, ["S1", "S2"])
        db.import_profile_bundle(conn, bundle(tmp_path, "S1"))
        custom = replace(ProfileConfig(), null_model=NullModelConfig(error_rate=0.02))
        with pytest.raises(ValueError, match="null_model.error_rate"):
            db.import_profile_bundle(conn, bundle(tmp_path, "S2", config=custom))
        assert conn.execute("SELECT count(*) FROM samples").fetchone()[0] == 1


@pytest.mark.parametrize("command", ["sync-profile", "profile-sra"])
def test_null_model_cli_overrides_one_toml_field_only(tmp_path, monkeypatch, command):
    path = tmp_path / "workflow.toml"
    path.write_text('[profile.null_model]\nerror_rate=0.02\np_threshold=0.001\nmax_total_reads=500\n')
    project = tmp_path / "project.duckdb"
    with db.connect(project) as conn:
        db.add_runs(conn, ["SRR1"])
    observed = []
    def fake_workflow(**kwargs):
        observed.append(kwargs["workflow_config"].profile.null_model)
        return workflows.SyncSummary(1, 0, 0, 0, 0)
    monkeypatch.setattr(workflows, "sync_remaining_profiles", fake_workflow)
    monkeypatch.setattr(workflows, "profile_sra_runs", fake_workflow)
    args = [command, "--db", str(project), "--cache-dir", str(tmp_path / "cache"),
            "--scratch-dir", str(tmp_path / "scratch"), "--output-dir", str(tmp_path / "outputs"),
            "--workflow-config", str(path), "--null-model-error-rate", "0.03"]
    if command == "profile-sra":
        remaining = tmp_path / "remaining.csv"
        remaining.write_text("run_id\nSRR1\n")
        args.extend(["--remaining-csv", str(remaining)])
    result = CliRunner().invoke(cli.cli, args)
    assert result.exit_code == 0, result.output
    assert observed == [NullModelConfig(error_rate=0.03, p_threshold=0.001, max_total_reads=500)]


def test_null_model_fallback_updates_asset_hash_and_checks_reuse(tmp_path):
    work = tmp_path / "work"; work.mkdir()
    (work / "profiling_contract.json").write_text('{"reference_hash":"ref"}')
    config = replace(ProfileConfig(), null_model=NullModelConfig(error_rate=0.02))
    calls = []
    class Runtime:
        def run(self, stage, command, **kwargs):
            calls.append(command)
            pl.DataFrame({"cov": [1], "max_error_count": [0]}).write_parquet(work / "null_model.parquet")
    kwargs = dict(profile_work=work, run_id="SRR1", logger=WorkflowLogger(), runtime=Runtime(), profile_config=config)
    model = workflows._ensure_null_model(**kwargs)
    assert calls[0][calls[0].index("--error-rate") + 1] == "0.02"
    assert json.loads((work / "profiling_contract.json").read_text())["null_model_hash"] == provenance.sha256_file(model)
    workflows._ensure_null_model(**kwargs)
    assert len(calls) == 1
    with pytest.raises(ValueError, match="cached null-model settings differ"):
        workflows._ensure_null_model(**{**kwargs, "profile_config": ProfileConfig()})


def test_custom_settings_never_relabel_unverified_cached_model(tmp_path):
    (tmp_path / "null_model.parquet").write_bytes(b"existing-model")
    with pytest.raises(ValueError, match="cannot verify custom null-model"):
        workflows._ensure_null_model(
            profile_work=tmp_path, run_id="SRR1", logger=WorkflowLogger(), runtime=None,
            profile_config=replace(ProfileConfig(), null_model=NullModelConfig(error_rate=0.02)),
        )
    assert (tmp_path / "null_model.parquet").read_bytes() == b"existing-model"


def test_explicit_defaults_reuse_the_same_default_model(tmp_path):
    model = tmp_path / "null_model.parquet"
    model.write_bytes(b"existing-default-model")
    config = replace(ProfileConfig(), null_model=NullModelConfig(**provenance.zipstrain_defaults()["null_model"]))
    assert workflows._ensure_null_model(
        profile_work=tmp_path, run_id="SRR1", logger=WorkflowLogger(), runtime=None,
        profile_config=config,
    ) == model
    assert model.read_bytes() == b"existing-default-model"


def test_real_zipstrain_generates_the_requested_null_model(tmp_path, monkeypatch):
    from polars.testing import assert_frame_equal
    from zipstrain.profile import _build_null_model_frame

    monkeypatch.setenv("PATH", str(Path(sys.executable).parent) + os.pathsep + os.environ["PATH"])
    config = replace(ProfileConfig(), null_model=NullModelConfig(error_rate=0.02, max_total_reads=20, p_threshold=0.001))
    result = workflows._ensure_null_model(
        profile_work=tmp_path, run_id="SRR1", logger=WorkflowLogger(),
        runtime=workflows._DirectRuntime(1), profile_config=config,
    )
    expected = _build_null_model_frame(0.02, 20, 0.001, "poisson")
    assert_frame_equal(pl.read_parquet(result), expected)
    assert json.loads((tmp_path / "null_model_settings.json").read_text()) == provenance.null_model_settings(config)


def test_fixed_reference_cache_and_commands_include_null_settings(tmp_path):
    reference, genes, stb = tmp_path / "reference.fa", tmp_path / "genes.fa", tmp_path / "reference.stb"
    reference.write_text(">chr\nACGT\n")
    genes.write_text(">chr_1 # 1 # 4 # 1\nACGT\n")
    stb.write_text("chr\tGCF_1.1\n")
    calls = []
    class Runtime:
        def threads(self, stage):
            return 1
        def run(self, stage, command, **kwargs):
            calls.append(command)
            if command[:3] == ["zipstrain", "utilities", "prepare_profiling"]:
                directory = Path(command[command.index("--output-dir") + 1]); directory.mkdir()
                for name in ("reference.fasta", "genomes_bed_file.bed", "gene_range_table.tsv", "null_model.parquet"):
                    (directory / name).write_text("fixture")
                (directory / "profiling_contract.json").write_text("{}")
            elif command[:2] == ["samtools", "faidx"]:
                Path(command[2] + ".fai").write_text("index")
            elif command[0] == "bowtie2-build":
                for path in workflows._bowtie_index_paths(Path(command[-1]))[:6]:
                    path.write_text("index")
            else:
                raise AssertionError(command)
    kwargs = dict(reference_genome=reference, reference_genome_genes=genes, reference_stb=stb,
                  cache_dir=tmp_path / "cache", runtime=Runtime(), logger=WorkflowLogger())
    default = workflows._prepare_fixed_reference_context(**kwargs)
    custom_config = replace(ProfileConfig(), null_model=NullModelConfig(error_rate=0.02))
    custom = workflows._prepare_fixed_reference_context(**kwargs, profile_config=custom_config)
    assert custom.reference_key != default.reference_key
    prepared = [command for command in calls if command[:3] == ["zipstrain", "utilities", "prepare_profiling"]]
    assert prepared[1][prepared[1].index("--error-rate") + 1] == "0.02"
    assert workflows._prepare_fixed_reference_context(**kwargs, profile_config=custom_config) == custom
    assert len([command for command in calls if command[:3] == ["zipstrain", "utilities", "prepare_profiling"]]) == 2


@pytest.mark.parametrize("value", ['null_model = "bad"', '[profile.null_model]\nunknown = 1'])
def test_bad_null_model_tables_rejected(tmp_path, value):
    path = tmp_path / "config.toml"
    path.write_text("[profile]\n" + value if not value.startswith("[profile.") else value)
    with pytest.raises(ValueError, match="null_model"):
        load_workflow_config(path, threads=1, sample_count=1)


@pytest.mark.parametrize("rows", [0, 1])
def test_missing_stats_diagnosed_without_partial_publication(tmp_path, rows):
    work = tmp_path / "scratch"; work.mkdir()
    output = tmp_path / "outputs"; output.mkdir()
    pl.DataFrame({"pos": list(range(rows))}, schema={"pos": pl.Int64}).write_parquet(work / "SRR1_profile.parquet")
    pl.DataFrame({"gene": ["g"]}).write_parquet(work / "SRR1_gene_stats.parquet")
    (work / "SRR1.bam").write_bytes(b"retained")
    expected = "zero rows" if rows == 0 else "read-location"
    with pytest.raises(FileNotFoundError, match=expected):
        workflows._publish_profile_outputs(run_id="SRR1", profile_work=work, output_dir=output)
    assert not list(output.iterdir())
    assert (work / "SRR1_profile.parquet").is_file()
    assert (work / "SRR1.bam").read_bytes() == b"retained"


def test_empty_output_does_not_publish_other_outputs(tmp_path):
    work = tmp_path / "work"; work.mkdir()
    output = tmp_path / "out"; output.mkdir()
    for name in ("profile", "gene_stats", "genome_stats"):
        (work / f"SRR1_{name}.parquet").write_bytes(b"" if name == "gene_stats" else b"data")
    with pytest.raises(FileNotFoundError, match="No bundle was published"):
        workflows._publish_profile_outputs(run_id="SRR1", profile_work=work, output_dir=output)
    assert not list(output.iterdir())


def test_partial_fastq_conversion_is_retried_not_reused(tmp_path, monkeypatch):
    scratch = tmp_path / "scratch"
    sample = scratch / "SRR1"
    archive = sample / "sra" / "SRR1" / "SRR1.sra"
    conversions = []
    def fake_run(command, **kwargs):
        if command[0] == "prefetch":
            archive.parent.mkdir(parents=True); archive.write_text("SRA")
        elif command[0] == "fasterq-dump":
            assert "--force" in command
            conversions.append(command)
            out = Path(command[command.index("--outdir") + 1])
            (out / "SRR1.fastq").write_text("partial")
            if len(conversions) == 1:
                raise subprocess.CalledProcessError(3, command, stderr="interrupted")
            (out / "SRR1.fastq").write_text("@r\nAC\n+\n!!\n")
        return subprocess.CompletedProcess(command, 0)
    monkeypatch.setattr(workflows.subprocess, "run", fake_run)
    kwargs = dict(run_ids=["SRR1"], db_file=tmp_path / "project.duckdb", cache_dir=tmp_path / "cache",
                  scratch_dir=scratch, raise_on_error=False, logger=WorkflowLogger())
    workflows.profile_sra_runs(**kwargs)
    assert archive.is_file()
    assert workflows._sample_fastqs(sample, required=False) == []
    assert (sample / "fastq_conversion.pending").exists()
    workflows.profile_sra_runs(**kwargs)
    assert len(conversions) == 2
    assert workflows._sample_fastqs(sample) == [sample / "reads" / "SRR1.fastq"]
    assert not (sample / "fastq_conversion.pending").exists()
    assert not archive.exists()

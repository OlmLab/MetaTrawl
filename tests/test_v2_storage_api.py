from dataclasses import replace
from pathlib import Path
import json
import shutil
import subprocess

import duckdb
import polars as pl
from polars.testing import assert_frame_equal
import pytest

from metatrawl import db, open_database, provenance, workflows, cache
from metatrawl.config import ProfileConfig, load_workflow_config
from zipstrain import gene_info as gi
from zipstrain.compare import PairComparison
from zipstrain.profile import prepare_profiling_assets
from zipstrain.profile_backends import ProfileRequest, run_profile


@pytest.fixture
def v2_inputs(tmp_path):
    reference = tmp_path / "ref.fna"
    reference.write_text(">chr1\nAAAGGATTTCCCAAAGGTTAA\n")
    genes = tmp_path / "genes.fna"
    genes.write_text(">chr1_1 # 1 # 9 # 1 # ID=1\nAAAGGATTT\n"
                     ">chr1_2 # 7 # 18 # -1 # ID=2\nAACCTTTGGGAA\n")
    stb = tmp_path / "ref.stb"
    stb.write_text("chr1\tGCF_1.1\n")
    assets = tmp_path / "assets"
    prepared = prepare_profiling_assets(reference_fasta=reference, gene_fasta=genes,
                                       stb_file=stb, output_dir=assets, max_total_reads=100)
    genome_cache = tmp_path / "cache"
    (genome_cache / "genomes").mkdir(parents=True)
    shutil.copy(reference, genome_cache / "genomes" / "GCF_1.1.fna")
    return reference, genes, stb, assets, prepared, genome_cache


def make_bundle(root, inputs, sample, *, collect=True, variant=True, backend="python"):
    reference, genes, stb, assets, prepared, _ = inputs
    sam, bam = root / f"{sample}.sam", root / f"{sample}.bam"
    with sam.open("w") as stream:
        stream.write("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:21\n")
        for i in range(15):
            seq = "AAAGGTTTTCCCAAAGGTTAA" if variant and i >= 10 else "AAAGGATTTCCCAAAGGTTAA"
            stream.write(f"r{i}\t0\tchr1\t1\t60\t21M\t*\t0\t0\t{seq}\t{'I'*21}\tNM:i:{int(variant and i>=10)}\n")
    subprocess.run(["samtools", "view", "-b", "-o", str(bam), str(sam)], check=True, capture_output=True)
    out = root / sample
    run_profile(ProfileRequest(
        bam_file=bam, bed_file=prepared.bed_file, stb_file=stb, null_model=assets / "null_model.parquet",
        reference_fasta=reference, gene_range_table=assets / "gene_info.parquet",
        profiling_contract=prepared.profiling_contract_file, output_dir=out,
        collect_codons=collect, allele_integration="weighted", num_chunks=1, max_concurrency=1,
        min_read_ani=0, read_inclusion="all-mapped",
    ), backend)
    profile = out / f"{sample}_profile.parquet"
    config = ProfileConfig(collect_codons=collect, allele_integration="weighted",
                           min_read_ani=0, read_inclusion="all-mapped", backend=backend)
    provenance.write_manifest(profile_file=profile, sample_id=sample, contract=provenance.profiling_contract(config),
                              references=provenance.reference_catalog(cache.PreparedReference(reference, genes, stb)),
                              details={"fixture": True})
    return db.ProfileBundle(
        sample, profile, out / f"{sample}_genome_stats.parquet",
        gene_stats_file=out / f"{sample}_gene_stats.parquet",
        gene_info_file=assets / "gene_info.parquet", reference_codons_file=assets / "reference_codons.parquet",
        codon_counts_file=out / f"{sample}_codon_counts.parquet" if collect else None,
    )


@pytest.mark.parametrize("mode", ["full", "allele-mask"])
def test_real_v2_bundle_retains_catalog_codons_and_selection(tmp_path, v2_inputs, mode):
    bundle = make_bundle(tmp_path, v2_inputs, "S1")
    path = tmp_path / "db.duckdb"
    expected = pl.read_parquet(bundle.gene_stats_file)
    counts = pl.read_parquet(bundle.codon_counts_file)
    with db.connect(path) as conn:
        db.configure_profile_storage(conn, mode=mode, min_cov=5 if mode == "allele-mask" else None)
        db.import_profile_bundle(conn, bundle, add_run_if_missing=True, cache_dir=v2_inputs[-1])
        assert conn.execute("SELECT schema_version FROM database_metadata").fetchone() == (2,)
        got = conn.execute('SELECT gene_id, "dN", "dS", "dN_dS", "pN", "pS", "pN_pS" FROM gene_stats ORDER BY gene_id').pl()
        assert_frame_equal(got, expected.select(got.columns).sort("gene_id"))
        assert got["gene_id"].dtype == pl.UInt64
        assert conn.execute("SELECT bool_and(collect_codons) FROM sample_profile_metadata").fetchone() == (True,)
        assert db.remaining_runs(conn) == []
    database = open_database(path)
    assert_frame_equal(database.codon_counts(sample_id="S1", genome="GCF_1.1").collect().drop("sample_id").sort(counts.columns),
                       counts.sort(counts.columns))
    assert database.genes(genome="GCF_1.1").collect().height == 2
    with db.connect(path) as conn:
        db.import_profile_bundle(conn, bundle, cache_dir=v2_inputs[-1])
        assert conn.execute("SELECT count(*) FROM gene_catalog").fetchone() == (2,)
        assert conn.execute("SELECT count(*) FROM codon_counts").fetchone() == (counts.height,)


def test_no_codons_keeps_selection_null_and_retains_empty_catalog_evidence(tmp_path, v2_inputs):
    bundle = make_bundle(tmp_path, v2_inputs, "S1", collect=False)
    with db.connect(tmp_path / "db.duckdb") as conn:
        db.import_profile_bundle(conn, bundle, add_run_if_missing=True)
        assert conn.execute('SELECT count(*) FROM gene_stats WHERE "dN" IS NOT NULL OR "pS" IS NOT NULL').fetchone() == (0,)
        assert conn.execute("SELECT count(*) FROM codon_counts").fetchone() == (0,)


def test_publish_discover_cleanup_and_partial_bundle_protection(tmp_path, v2_inputs):
    bundle = make_bundle(tmp_path, v2_inputs, "S1")
    outputs = tmp_path / "published"
    outputs.mkdir()
    workflows._publish_profile_outputs(run_id="S1", profile_work=bundle.profile_file.parent,
                                      output_dir=outputs, gene_info_file=bundle.gene_info_file, collect_codons=True)
    provenance.write_manifest(profile_file=outputs/"S1.profile.parquet", sample_id="S1",
                              contract=provenance.profiling_contract(ProfileConfig()),
                              references={}, details={})
    assert workflows._published_profile_bundle_complete(run_id="S1", output_dir=outputs, require_sylph=False)
    imported = workflows.discover_profile_bundle(output_dir=outputs, run_id="S1")
    assert imported.codon_counts_file.is_file()
    imported.codon_counts_file.unlink()
    assert not workflows._published_profile_bundle_complete(run_id="S1", output_dir=outputs, require_sylph=False)
    with db.connect(tmp_path / "db.duckdb") as conn:
        with pytest.raises(ValueError, match="retained codon-count"):
            db.import_profile_bundle(conn, replace(imported, codon_counts_file=None), add_run_if_missing=True)
        assert conn.execute("SELECT count(*) FROM samples").fetchone() == (0,)
    workflows.cleanup_profile_bundle(imported)
    assert not list(outputs.glob("*.parquet"))
    assert bundle.gene_info_file.is_file()
    assert bundle.reference_codons_file.is_file()  # Shared assets survive publication and cleanup.


@pytest.mark.parametrize("mode", ["full", "allele-mask"])
def test_pair_api_matches_zipstrain_and_exports(tmp_path, v2_inputs, mode):
    a = make_bundle(tmp_path, v2_inputs, "S1", variant=False)
    b = make_bundle(tmp_path, v2_inputs, "S2")
    path = tmp_path / "db.duckdb"
    with db.connect(path) as conn:
        db.configure_profile_storage(conn, mode=mode, min_cov=5 if mode == "allele-mask" else None)
        db.import_profile_bundles(conn, [a, b], add_runs_if_missing=True, cache_dir=v2_inputs[-1])
    pair = open_database(path).pair("S1", "S2", genome="GCF_1.1")
    expected = PairComparison(a.profile_file, b.profile_file, sample_1="S1", sample_2="S2")
    assert_frame_equal(pair.compare().collect(), expected.compare("genome_ani+max_identical_stretch").collect())
    assert_frame_equal(pair.identical_stretches(min_length=1).collect(), expected.identical_stretches(min_length=1).collect())
    assert_frame_equal(pair.mosaic_regions().collect(), expected.mosaic_regions().collect())
    query = pair.identical_stretches(min_length=1)
    output = tmp_path / "stretches.parquet"
    query.sink_parquet(output)
    assert_frame_equal(pl.read_parquet(output), query.collect())
    with pytest.raises(FileExistsError):
        query.sink_parquet(output)
    if mode == "allele-mask":
        with pytest.raises(ValueError, match="fixed min_cov"):
            open_database(path).pair("S1", "S2", genome="GCF_1.1", min_cov=6).compare().collect()


def test_v1_database_is_readable_but_not_writable_or_relabelled(tmp_path):
    path = tmp_path / "v1.duckdb"
    with db.connect(path) as conn:
        conn.execute("UPDATE database_metadata SET schema_version=1")
    assert open_database(path).metadata()["schema_version"] == 1
    with pytest.raises(ValueError, match="schema"):
        db.connect(path)
    with duckdb.connect(str(path), read_only=True) as conn:
        assert conn.execute("SELECT schema_version FROM database_metadata").fetchone() == (1,)


def test_tampered_codons_roll_back_entire_import(tmp_path, v2_inputs):
    bundle = make_bundle(tmp_path, v2_inputs, "S1")
    source = bundle.codon_counts_file
    meta = pl.read_parquet_metadata(source)
    meta = {k: v for k,v in meta.items() if k.startswith("zipstrain_")}
    meta["zipstrain_source_profile_sha256"] = "wrong"
    pl.read_parquet(source).write_parquet(source, metadata=meta)
    with db.connect(tmp_path / "db.duckdb") as conn:
        with pytest.raises(ValueError, match="codon evidence does not match"):
            db.import_profile_bundle(conn, bundle, add_run_if_missing=True)
        assert conn.execute("SELECT count(*) FROM gene_catalog").fetchone() == (0,)
        assert conn.execute("SELECT count(*) FROM codon_counts").fetchone() == (0,)


def test_selection_config_validation_and_defaults(tmp_path):
    path = tmp_path / "config.toml"
    path.write_text('[profile]\ncollect_codons=true\nallele_integration="weighted"\n')
    config = load_workflow_config(path, threads=1, sample_count=1).profile
    assert config.collect_codons and config.allele_integration == "weighted"
    assert config.min_freq == .01
    with pytest.raises(ValueError, match="collect_codons"):
        ProfileConfig(collect_codons="false")
    with pytest.raises(ValueError, match="allele_integration"):
        ProfileConfig(allele_integration="bad")


@pytest.mark.parametrize("command", ["sync-profile", "profile-sra"])
def test_selection_cli_overrides_toml(tmp_path, monkeypatch, command):
    from click.testing import CliRunner
    from metatrawl.cli import cli
    config = tmp_path / "workflow.toml"
    config.write_text('[profile]\ncollect_codons=true\nallele_integration="consensus"\n')
    database = tmp_path / "project.duckdb"
    with db.connect(database) as conn:
        db.add_runs(conn, ["S1"])
    captured = []
    def run(**kwargs):
        captured.append(kwargs["workflow_config"].profile)
        return workflows.SyncSummary(1, 0, 0, 0, 0)
    monkeypatch.setattr(workflows, "sync_remaining_profiles", run)
    monkeypatch.setattr(workflows, "profile_sra_runs", run)
    args = [command, "--db", str(database), "--cache-dir", str(tmp_path/"cache"),
            "--scratch-dir", str(tmp_path/"scratch"), "--output-dir", str(tmp_path/"outputs"),
            "--workflow-config", str(config), "--no-collect-codons", "--allele-integration", "weighted"]
    if command == "profile-sra":
        pending = tmp_path/"pending.csv"
        pending.write_text("run_id\nS1\n")
        args += ["--remaining-csv", str(pending)]
    result = CliRunner().invoke(cli, args)
    assert result.exit_code == 0, result.output
    assert not captured[0].collect_codons and captured[0].allele_integration == "weighted"


@pytest.mark.parametrize("backend", ["python", "rust_profiler"])
def test_real_worker_publishes_evidence_before_scratch_cleanup(tmp_path, v2_inputs, backend):
    from click.testing import CliRunner
    from zipstrain.cli import cli as zip_cli
    from metatrawl.logging import WorkflowLogger
    reference, genes, stb, assets, prepared, _ = v2_inputs
    seed = make_bundle(tmp_path, v2_inputs, "S1", backend=backend)
    scratch = tmp_path / "scratch" / "S1"
    scratch.mkdir(parents=True)
    shutil.copy(tmp_path/"S1.bam", scratch/"S1.bam")
    (scratch/"read_layout.txt").write_text("single\n")
    config = ProfileConfig(collect_codons=True, allele_integration="weighted", backend=backend,
                           min_read_ani=0, read_inclusion="all-mapped")
    (scratch/"profiling_settings.json").write_text(json.dumps(provenance.profiling_contract(config)))
    commands = []
    class Runtime:
        def threads(self, stage):
            return 1
        def run(self, stage, command, **kwargs):
            commands.append(command)
            if command[0] == "zipstrain":
                result = CliRunner().invoke(zip_cli, command[1:])
                if result.exit_code:
                    raise RuntimeError(result.output + repr(result.exception))
            else:
                subprocess.run(command, check=True, capture_output=True)
    outputs = tmp_path/"outputs"
    outputs.mkdir()
    (outputs/"S1.sylph.tsv").write_text("Genome_file\tTaxonomic_abundance\nGCF_1.1\t100\n")
    workflows._run_alignment_and_profile(
        run_id="S1", sample_scratch=scratch, reference=cache.PreparedReference(reference, genes, stb),
        output_dir=outputs, runtime=Runtime(), profile_config=config, logger=WorkflowLogger())
    command = next(command for command in commands if "profile-single" in command)
    assert "--collect-codons" in command
    assert command[command.index("--allele-integration")+1] == "weighted"
    assert not (scratch/"S1.bam").exists()
    assert not (scratch/"zipstrain_profile").exists()
    bundle = workflows.discover_profile_bundle(output_dir=outputs, run_id="S1")
    assert bundle.codon_counts_file.is_file()
    with db.connect(tmp_path/"worker.duckdb") as conn:
        db.import_profile_bundle(conn, bundle, add_run_if_missing=True)
        assert conn.execute("SELECT count(*) FROM codon_counts").fetchone()[0] > 0
    workflows.cleanup_profile_bundle(bundle)
    assert not list(outputs.glob("*.parquet"))


def test_pair_rejects_incompatible_metadata_before_reading_positions(tmp_path, v2_inputs):
    a = make_bundle(tmp_path, v2_inputs, "S1")
    b = make_bundle(tmp_path, v2_inputs, "S2")
    path = tmp_path/"db.duckdb"
    with db.connect(path) as conn:
        db.import_profile_bundles(conn, [a,b], add_runs_if_missing=True)
        row = conn.execute("SELECT manifest_json FROM sample_provenance WHERE sample_id='S2'").fetchone()
        manifest = json.loads(row[0])
        manifest["references"]["GCF_1.1"]["sequence_hash"] = "different"
        conn.execute("UPDATE sample_provenance SET manifest_json=? WHERE sample_id='S2'", [json.dumps(manifest)])
        conn.execute("DROP TABLE profile_positions")
    with pytest.raises(ValueError, match="genome reference"):
        open_database(path).pair("S1", "S2", genome="GCF_1.1").compare().collect()


@pytest.mark.parametrize("mode", ["full", "allele-mask"])
@pytest.mark.parametrize("integration", ["consensus", "weighted"])
def test_pair_selection_matches_zipstrain_without_position_access(tmp_path, v2_inputs, mode, integration):
    from zipstrain import codons
    a = make_bundle(tmp_path, v2_inputs, "S1", variant=False)
    b = make_bundle(tmp_path, v2_inputs, "S2")
    expected = codons.pairwise_selection(a.codon_counts_file, b.codon_counts_file, a.gene_info_file,
                                        allele_integration=integration).collect().sort("gene_id")
    path = tmp_path/"project.duckdb"
    with db.connect(path) as conn:
        db.configure_profile_storage(conn, mode=mode, min_cov=5 if mode == "allele-mask" else None)
        db.import_profile_bundles(conn, [a, b], add_runs_if_missing=True, cache_dir=v2_inputs[-1])
        # Selection must use read-phased codons, never reconstructed position counts.
        conn.execute("DROP TABLE profile_positions")
    pair = open_database(path).pair("S1", "S2", genome="GCF_1.1")
    query = pair.gene_selection(allele_integration=integration)
    result = query.collect().sort("gene_id")
    assert_frame_equal(result.select(expected.columns), expected)
    assert result["gene"].null_count() == 0
    assert result["sample_1"].to_list() == ["S1"] * expected.height
    output = tmp_path/"selection.parquet"
    query.sink_parquet(output)
    assert_frame_equal(pl.read_parquet(output).sort("gene_id"), result)
    assert_frame_equal(pair.compare("gene_selection", allele_integration=integration).collect().sort("gene_id"), result)


def test_pair_selection_requires_retained_codons(tmp_path, v2_inputs):
    a = make_bundle(tmp_path, v2_inputs, "S1", collect=False)
    b = make_bundle(tmp_path, v2_inputs, "S2", collect=False)
    path = tmp_path/"project.duckdb"
    with db.connect(path) as conn:
        db.import_profile_bundles(conn, [a, b], add_runs_if_missing=True)
    with pytest.raises(ValueError, match="collect_codons=true"):
        open_database(path).pair("S1", "S2", genome="GCF_1.1").gene_selection().collect()


def test_matrix_compare_does_not_accept_pairwise_selection():
    with pytest.raises(ValueError, match="API"):
        workflows.normalize_calculations("gene_selection")


def test_pair_selection_with_different_sample_reference_subsets(tmp_path, v2_inputs):
    from zipstrain import codons
    a = make_bundle(tmp_path, v2_inputs, "S1", variant=False)
    baseline = make_bundle(tmp_path, v2_inputs, "baseline")
    expected = codons.pairwise_selection(a.codon_counts_file, baseline.codon_counts_file,
                                        a.gene_info_file, allele_integration="weighted").collect()
    root = tmp_path/"expanded"
    root.mkdir()
    reference, genes, stb = root/"ref.fna", root/"genes.fna", root/"ref.stb"
    reference.write_text(v2_inputs[0].read_text() + ">chr2\nAAAGGATTTCCCAAAGGTTAA\n")
    genes.write_text(v2_inputs[1].read_text() + ">chr2_1 # 1 # 9 # 1 # ID=3\nAAAGGATTT\n")
    stb.write_text(v2_inputs[2].read_text() + "chr2\tGCF_2.1\n")
    assets = root/"assets"
    prepared = prepare_profiling_assets(reference_fasta=reference, gene_fasta=genes,
                                       stb_file=stb, output_dir=assets, max_total_reads=100)
    b = make_bundle(root, (reference, genes, stb, assets, prepared, v2_inputs[-1]), "S2")
    path = tmp_path/"project.duckdb"
    with db.connect(path) as conn:
        db.import_profile_bundles(conn, [a, b], add_runs_if_missing=True)
        assert conn.execute("SELECT count(*) FROM gene_catalog").fetchone() == (3,)
    result = open_database(path).pair("S1", "S2", genome="GCF_1.1").gene_selection(
        allele_integration="weighted").collect()
    assert_frame_equal(result.select(expected.columns).sort("gene_id"), expected.sort("gene_id"))


def test_storage_conversion_preserves_selection_evidence(tmp_path, v2_inputs):
    from metatrawl import migration
    a = make_bundle(tmp_path, v2_inputs, "S1", variant=False)
    b = make_bundle(tmp_path, v2_inputs, "S2")
    source, target = tmp_path/"full.duckdb", tmp_path/"mask.duckdb"
    with db.connect(source) as conn:
        db.import_profile_bundles(conn, [a, b], add_runs_if_missing=True)
    summary = migration.migrate_full_database(source_db=source, output_db=target,
                                             cache_dir=v2_inputs[-1], min_cov=5)
    assert summary.completed_samples == 2
    expected = open_database(source).pair("S1", "S2", genome="GCF_1.1").gene_selection().collect()
    got = open_database(target).pair("S1", "S2", genome="GCF_1.1").gene_selection().collect()
    assert_frame_equal(got.sort("gene_id"), expected.sort("gene_id"))


@pytest.mark.parametrize("sparse", [False, True])
def test_v2_gene_catalog_matrix_build_and_comparison(tmp_path, v2_inputs, sparse):
    from metatrawl import matrix_hdf5
    from zipstrain import matrix_pairs
    a = make_bundle(tmp_path, v2_inputs, "S1", variant=False)
    b = make_bundle(tmp_path, v2_inputs, "S2")
    project = tmp_path/"project.duckdb"
    with db.connect(project) as conn:
        db.import_profile_bundles(conn, [a, b], add_runs_if_missing=True)
        matrix_hdf5.build_matrix_hdf5_from_duckdb(
            conn, sample_ids=["S1", "S2"], output_file=tmp_path/"matrix.h5", genome="GCF_1.1",
            bed_file=v2_inputs[4].bed_file, stb_file=v2_inputs[2], gene_range_table=a.gene_info_file,
            sparse=sparse)
    matrix_pairs.matrix_compare(matrix_db_file=tmp_path/"matrix.h5", output_file=tmp_path/"compare.duckdb",
                                calculate="all", backend="numpy")
    with duckdb.connect(str(tmp_path/"compare.duckdb"), read_only=True) as conn:
        genome = conn.execute("SELECT genome_ani, max_identical_stretch FROM matrix_compare_results").pl()
        expected = PairComparison(a.profile_file, b.profile_file).compare("genome_ani+max_identical_stretch").collect()
        assert_frame_equal(genome, expected.select(genome.columns))
        assert conn.execute("SELECT count(*) FROM matrix_compare_gene_results").fetchone() == (2,)
        columns = {row[1] for row in conn.execute("PRAGMA table_info('matrix_compare_gene_results')").fetchall()}
        assert not columns & {"dN", "dS", "pN", "pS"}


def test_invalid_gene_stats_roll_back(tmp_path, v2_inputs):
    bundle = make_bundle(tmp_path, v2_inputs, "S1")
    source = bundle.gene_stats_file
    from metatrawl.v2_storage import metadata
    pl.read_parquet(source).with_columns(start=pl.col("start") + 1).write_parquet(source, metadata=metadata(source))
    with db.connect(tmp_path/"project.duckdb") as conn:
        with pytest.raises(ValueError, match="gene statistics do not match"):
            db.import_profile_bundle(conn, bundle, add_run_if_missing=True)
        assert conn.execute("SELECT count(*) FROM gene_catalog").fetchone() == (0,)
        assert conn.execute("SELECT count(*) FROM samples").fetchone() == (0,)

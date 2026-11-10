from pathlib import Path

import duckdb
import polars as pl
from polars.testing import assert_frame_equal
import pytest

from metatrawl import allele_mask, open_database
from metatrawl import db as registry
from metatrawl.api import ProfileReferenceUnavailableError


BASES = "ACGT"
KEY = ["sample_id", "genome", "chrom", "pos"]
MIN_COV = 7


@pytest.fixture
def snp_databases(tmp_path):
    paths = {mode: tmp_path / f"{mode}.duckdb" for mode in ("full", "allele-mask")}
    cache = tmp_path / "cache"
    (cache / "genomes").mkdir(parents=True)
    # All nonempty allele sets against every reference base, plus ambiguous,
    # uncovered, and subthreshold positions. Use an odd-sized scaffold.
    sequence = "ACGTN" * 15 + "ACA"
    for genome in ("genome_1", "genome_2"):
        (cache / "genomes" / f"{genome}.fna").write_text(
            f">contig_1\n{sequence}\n>contig_2\n{sequence}\n"
        )
    expected = []
    for sample in ("sample_a", "sample'b"):
        rows = []
        for genome in ("genome_1", "genome_2"):
            for chrom in ("contig_1", "contig_2"):
                for pos, ref in enumerate(sequence, 1):
                    mask = ((pos - 1) // 5 + (sample == "sample'b")) % 15 + 1
                    counts = [8 + i if mask & (1 << i) else 0 for i in range(4)]
                    if pos == len(sequence) - 2:
                        counts = [0, 0, 0, 0]
                    elif pos == len(sequence) - 1:
                        counts = [1, 1, 1, 1]
                    elif pos == len(sequence):
                        counts = [0, 0, 0, 7]
                    ref_bit = 1 << BASES.index(ref) if ref in BASES else 0
                    rows.append((chrom, genome, pos, *counts, ref_bit))
                    observed = [base for base, count in zip(BASES, counts) if count]
                    alts = [base for base in observed if base != ref]
                    if ref_bit and sum(counts) >= MIN_COV and alts:
                        expected.append((sample, genome, chrom, pos, ref, alts, len(observed), ref in observed))
        frame = pl.DataFrame(
            rows,
            schema=["chrom", "genome", "pos", *BASES, "ref_base_bitmask"],
            orient="row",
        )
        profile = tmp_path / f"{sample}.parquet"
        frame.write_parquet(profile)
        with registry.connect(paths["full"]) as conn:
            registry._insert_profile_positions(conn, sample_id=sample, profile_file=profile)
        with registry.connect(paths["allele-mask"]) as conn:
            registry.configure_profile_storage(conn, mode="allele-mask", min_cov=MIN_COV)
            allele_mask.store_profile_parquet(
                conn, sample_id=sample, profile_file=profile, min_cov=MIN_COV, cache_dir=cache,
            )
    # SNP queries must not need an external FASTA after import.
    for fasta in (cache / "genomes").glob("*.fna"):
        fasta.unlink()
    return {mode: open_database(path) for mode, path in paths.items()}, sorted(expected)


def test_reference_snp_sites_and_alleles_agree_across_storage_modes(snp_databases):
    dbs, expected = snp_databases
    compact = dbs["allele-mask"]
    full = dbs["full"]

    for sample in ("sample_a", "sample'b"):
        observed = compact.sample(sample).snps().collect().sort(KEY)
        assert observed.rows() == [row for row in expected if row[0] == sample]
        assert_frame_equal(
            observed,
            full.sample(sample).snps(min_cov=MIN_COV).collect().select(observed.columns).sort(KEY),
        )
    assert len(expected) > 300


@pytest.mark.parametrize("mode", ["full", "allele-mask"])
def test_snp_query_scope_lazy_and_streaming_parquet(snp_databases, tmp_path, mode):
    dbs, expected = snp_databases
    db = dbs[mode]
    query = db.genome("genome_2").snps(chrom="contig_2", min_cov=MIN_COV)
    result = query.collect().sort(KEY)
    assert result["genome"].unique().to_list() == ["genome_2"]
    assert result["chrom"].unique().to_list() == ["contig_2"]
    assert result.height == sum(row[1:3] == ("genome_2", "contig_2") for row in expected)
    assert_frame_equal(
        query.lazy().filter(pl.col("sample_id") == "sample'b").collect().sort(KEY),
        db.snps(sample_id="sample'b", genome="genome_2", chrom="contig_2", min_cov=MIN_COV).collect().sort(KEY),
    )

    file = tmp_path / "nested" / "snps.parquet"
    assert query.sink_parquet(file) == file
    assert_frame_equal(pl.read_parquet(file).sort(KEY), result)
    with pytest.raises(FileExistsError):
        query.sink_parquet(file)
    empty_query = db.sample("missing").snps()
    assert empty_query.collect().is_empty()
    empty_query.sink_parquet(file, overwrite=True)
    assert_frame_equal(pl.read_parquet(file), result.head(0))


def test_full_snp_details_match_zipstrain(snp_databases):
    from zipstrain.profile import get_reference_snps

    dbs, _ = snp_databases
    sample = dbs["full"].sample("sample_a")
    result = sample.snps("genome_1", min_cov=MIN_COV).collect()
    expected = get_reference_snps(
        sample.profile("genome_1").lazy().with_columns(
            pl.col(*BASES).cast(pl.Int64), pl.lit("NA").alias("gene"),
        ),
        min_cov=MIN_COV,
    ).collect()
    columns = ["chrom", "pos", "ref_base", "con_base", "class", "ref_freq", "con_freq", "position_coverage"]
    assert_frame_equal(
        result.select(columns).sort("chrom", "pos"),
        expected.select(columns).sort("chrom", "pos"),
    )
    assert set(result["class"]) == {"SNS", "SNV", "con_SNV", "pop_SNV"}
    assert (result["alt_freq"] + result["ref_freq"]).min() == pytest.approx(1.0)


def test_full_snp_tie_break_depth_and_missing_reference(tmp_path):
    path = tmp_path / "full.duckdb"
    with registry.connect(path) as conn:
        conn.execute("""
            INSERT INTO profile_positions VALUES
                ('s', 'c', 1, 'g', 65535, 65535, 0, 0, 2),
                ('s', 'c', 2, 'g', 5, 0, 0, 0, 2),
                ('s', 'c', 3, 'g', 5, 0, 0, 0, NULL),
                ('s', 'c', 4, 'g', 5, 0, 0, 0, 3),
                ('s', 'c', 5, 'g', 4, 0, 0, 0, 2)
        """)
    snps = open_database(path).sample("s").snps().collect().sort("pos")
    assert snps.select("pos", "con_base", "position_coverage", "class").rows() == [
        (1, "A", 131070, "con_SNV"),
        (2, "A", 5, "SNS"),
    ]


def test_snp_validation_and_fixed_coverage(snp_databases):
    dbs, _ = snp_databases
    compact = dbs["allele-mask"]
    with pytest.raises(ValueError, match="requires a sample_id or genome"):
        compact.snps()
    for value in (True, 0, -1, 1.5, "7"):
        with pytest.raises(ValueError, match="positive integer"):
            compact.sample("sample_a").snps(min_cov=value)
    for value in (MIN_COV - 1, MIN_COV + 1):
        with pytest.raises(ValueError, match=f"fixed min_cov={MIN_COV}"):
            compact.sample("sample_a").snps(min_cov=value)


def test_legacy_full_schema_is_read_without_migration(tmp_path):
    path = tmp_path / "legacy.duckdb"
    with duckdb.connect(str(path)) as conn:
        conn.execute("""
            CREATE TABLE profile_positions AS
            SELECT 's' AS sample_id, 'g' AS genome, 'c' AS chrom, 1::BIGINT AS pos,
                   0.0::DOUBLE AS A, 5.0::DOUBLE AS C, 0.0::DOUBLE AS G, 0.0::DOUBLE AS T
        """)
    with pytest.raises(ProfileReferenceUnavailableError, match="reference FASTA"):
        open_database(path).sample("s").snps()
    with duckdb.connect(str(path)) as conn:
        conn.execute("ALTER TABLE profile_positions ADD COLUMN ref_base_bitmask UTINYINT DEFAULT 1")
    assert open_database(path).sample("s").snps().collect()["alt_bases"].to_list() == [["C"]]
    with duckdb.connect(str(path), read_only=True) as conn:
        assert conn.execute("SHOW TABLES").fetchall() == [("profile_positions",)]


def test_compact_snp_export_does_not_publish_corrupt_results(snp_databases, tmp_path):
    dbs, _ = snp_databases
    compact = dbs["allele-mask"]
    with registry.connect(compact.path) as conn:
        conn.execute("UPDATE allele_mask_profile_blocks SET payload_hash = 'invalid' WHERE sample_id = 'sample_a'")
    file = tmp_path / "snps.parquet"
    file.write_bytes(b"keep the old output")
    with pytest.raises(ValueError, match="checksum mismatch for sample_a"):
        compact.sample("sample_a").snps().sink_parquet(file, overwrite=True)
    assert file.read_bytes() == b"keep the old output"
    assert not list(tmp_path.glob(".snps.parquet.*"))

from __future__ import annotations

from pathlib import Path

import polars as pl
import pytest

from metatrawl import MetaTrawlDatabase, open_database
from metatrawl import db as registry


def _database(tmp_path: Path) -> Path:
    db_path = tmp_path / "metatrawl.duckdb"
    with registry.connect(db_path) as conn:
        conn.executemany(
            "INSERT INTO samples VALUES (?, ?, 'complete', 1, 1)",
            [("sample_a", "run_a"), ("sample_b", "run_b")],
        )
        conn.executemany(
            "INSERT INTO genome_stats VALUES (?, ?, ?, ?, ?)",
            [
                ("sample_a", "genome_1", 5.0, 0.9, 0.8),
                ("sample_a", "genome_2", 2.0, 0.5, 0.4),
                ("sample_b", "genome_1", 7.0, 0.95, 0.85),
            ],
        )
        conn.executemany(
            "INSERT INTO gene_stats VALUES (?, ?, ?, ?, ?, ?)",
            [
                ("sample_a", "genome_1", "gene_1", 5.0, 0.9, 0.8),
                ("sample_a", "genome_1", "gene_2", 4.0, 0.8, 0.7),
                ("sample_b", "genome_1", "gene_1", 7.0, 0.95, 0.85),
            ],
        )
        conn.executemany(
            "INSERT INTO profile_positions VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
            [
                ("sample_a", "contig_1", 1, "genome_1", "gene_1", 5.0, 0.0, 0.0, 0.0),
                ("sample_a", "contig_1", 2, "genome_1", "gene_1", 0.0, 5.0, 0.0, 0.0),
                ("sample_b", "contig_1", 1, "genome_1", "gene_1", 6.0, 0.0, 0.0, 0.0),
            ],
        )
        conn.executemany(
            "INSERT INTO sylph_abundance VALUES (?, ?, ?, ?)",
            [
                ("sample_a", "genome_1", "GCF_1", 0.4),
                ("sample_b", "genome_1", "GCF_1", 0.6),
            ],
        )
    return db_path


def test_genome_view_queries_stats_across_samples(tmp_path: Path) -> None:
    database = open_database(_database(tmp_path))

    genome_stats = database.genome("genome_1").genome_stats().collect()
    gene_stats = database.genome("genome_1").gene_stats("gene_1").collect()

    assert genome_stats["sample_id"].to_list() == ["sample_a", "sample_b"]
    assert gene_stats.select("sample_id", "gene").rows() == [
        ("sample_a", "gene_1"),
        ("sample_b", "gene_1"),
    ]


def test_sample_view_queries_profiles_genomes_and_genes(tmp_path: Path) -> None:
    sample = MetaTrawlDatabase(_database(tmp_path)).sample("sample_a")

    assert sample.genome_stats().collect()["genome"].to_list() == ["genome_1", "genome_2"]
    assert sample.gene_stats(genome="genome_1").collect()["gene"].to_list() == ["gene_1", "gene_2"]
    profile = sample.profile(genome="genome_1", gene="gene_1").collect()
    assert profile.shape == (2, 9)
    assert profile["pos"].to_list() == [1, 2]
    assert profile.select("A", "C", "G", "T").schema == {
        "A": pl.UInt16,
        "C": pl.UInt16,
        "G": pl.UInt16,
        "T": pl.UInt16,
    }


def test_query_lazy_supports_polars_expressions(tmp_path: Path) -> None:
    query = open_database(_database(tmp_path)).genome("genome_1").genome_stats()

    result = query.lazy().filter(pl.col("coverage") >= 6).collect()

    assert result["sample_id"].to_list() == ["sample_b"]


def test_query_sinks_directly_to_parquet(tmp_path: Path) -> None:
    query = open_database(_database(tmp_path)).genome("genome_1").profiles()
    output = tmp_path / "exports" / "profiles.parquet"

    returned = query.sink_parquet(output)

    assert returned == output.resolve()
    exported = pl.read_parquet(output)
    assert exported.shape == (3, 9)
    assert set(exported["sample_id"]) == {"sample_a", "sample_b"}
    with pytest.raises(FileExistsError):
        query.sink_parquet(output)


def test_query_parquet_sink_can_overwrite(tmp_path: Path) -> None:
    query = open_database(_database(tmp_path)).sample("sample_a").genome_stats()
    output = tmp_path / "stats.parquet"
    output.write_text("old")

    query.sink_parquet(output, overwrite=True)

    assert pl.read_parquet(output)["genome"].to_list() == ["genome_1", "genome_2"]


def test_missing_database_and_empty_keys_are_rejected(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError):
        open_database(tmp_path / "missing.duckdb")

    database = open_database(_database(tmp_path))
    with pytest.raises(ValueError, match="genome cannot be empty"):
        database.genome(" ")
    with pytest.raises(ValueError, match="sample_id cannot be empty"):
        database.sample("")


def test_genome_view_queries_sylph_by_accession(tmp_path: Path) -> None:
    result = open_database(_database(tmp_path)).genome("GCF_1").sylph_abundance().collect()

    assert result["sample_id"].to_list() == ["sample_a", "sample_b"]


def test_connect_migrates_legacy_float_profile_counts(tmp_path: Path) -> None:
    db_path = tmp_path / "legacy.duckdb"
    with registry.connect(db_path) as conn:
        conn.execute(
            """
            CREATE OR REPLACE TABLE profile_positions (
                sample_id VARCHAR,
                chrom VARCHAR,
                pos BIGINT,
                genome VARCHAR,
                gene VARCHAR,
                A DOUBLE,
                C DOUBLE,
                G DOUBLE,
                T DOUBLE
            )
            """
        )
        conn.execute(
            "INSERT INTO profile_positions VALUES ('sample', 'contig', 1, 'genome', 'gene', 5.0, 0.0, 1.0, 2.0)"
        )

    with registry.connect(db_path) as conn:
        types = dict(
            conn.execute(
                """
                SELECT column_name, data_type
                FROM information_schema.columns
                WHERE table_name = 'profile_positions'
                  AND column_name IN ('A', 'C', 'G', 'T')
                """
            ).fetchall()
        )
        values = conn.execute("SELECT A, C, G, T FROM profile_positions").fetchone()

    assert types == {"A": "USMALLINT", "C": "USMALLINT", "G": "USMALLINT", "T": "USMALLINT"}
    assert values == (5, 0, 1, 2)


def test_connect_normalizes_existing_sylph_genome_paths(tmp_path: Path) -> None:
    db_path = _database(tmp_path)
    with registry.connect(db_path) as conn:
        conn.execute(
            """
            UPDATE sylph_abundance
            SET genome = ?, accession = ?
            WHERE sample_id = 'sample_a'
            """,
            [
                "gtdb/database/GCF/901/875/305/GCF_901875305.1_genomic.fna.gz",
                "gtdb/database/GCF/901/875/305/GCF_901875305.1_genomic.fna.gz",
            ],
        )

    with registry.connect(db_path) as conn:
        value = conn.execute(
            "SELECT genome, accession FROM sylph_abundance WHERE sample_id = 'sample_a'"
        ).fetchone()

    assert value == ("GCF_901875305.1", "GCF_901875305.1")


def test_database_lists_distinct_genomes_and_samples(tmp_path: Path) -> None:
    database = open_database(_database(tmp_path))

    assert database.genomes().collect()["genome"].to_list() == ["genome_1", "genome_2"]
    samples = database.samples().collect()
    assert samples["sample_id"].to_list() == ["sample_a", "sample_b"]
    assert samples["status"].to_list() == ["complete", "complete"]


def test_database_filters_genomes_and_samples_with_regex(tmp_path: Path) -> None:
    database = open_database(_database(tmp_path))

    assert database.genomes(pattern=r"_2$").collect()["genome"].to_list() == ["genome_2"]
    assert database.samples(pattern=r"sample_[ab]").collect()["sample_id"].to_list() == [
        "sample_a",
        "sample_b",
    ]

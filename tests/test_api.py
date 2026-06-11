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

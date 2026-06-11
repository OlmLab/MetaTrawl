"""Public Python API for querying a MetaTrawl project database."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import duckdb
import polars as pl


class MetaTrawlDatabase:
    """Read data from a MetaTrawl DuckDB database."""

    def __init__(self, db_path: str | Path) -> None:
        self.path = Path(db_path).expanduser().resolve()
        if not self.path.is_file():
            raise FileNotFoundError(f"MetaTrawl database does not exist: {self.path}")

    def genome(self, genome: str) -> GenomeView:
        """Return a genome-centric view across all samples."""
        return GenomeView(self.path, _required_text(genome, "genome"))

    def sample(self, sample_id: str) -> SampleView:
        """Return a view of all stored data for one sample."""
        return SampleView(self.path, _required_text(sample_id, "sample_id"))


@dataclass(frozen=True)
class Query:
    """A reusable, parameterized query against a MetaTrawl database."""

    db_path: Path
    sql: str
    parameters: tuple[Any, ...] = ()

    def collect(self) -> pl.DataFrame:
        """Execute the query and return an eager Polars DataFrame."""
        with duckdb.connect(str(self.db_path), read_only=True) as conn:
            return conn.execute(self.sql, list(self.parameters)).pl()

    def lazy(self) -> pl.LazyFrame:
        """Return the query result for additional lazy Polars operations.

        The DuckDB query is materialized when this method is called. Use
        :meth:`sink_parquet` for a streaming export that avoids Python memory.
        """
        return self.collect().lazy()

    def sink_parquet(
        self,
        output_file: str | Path,
        *,
        compression: str = "zstd",
        overwrite: bool = False,
    ) -> Path:
        """Stream the query result directly from DuckDB to Parquet."""
        path = Path(output_file).expanduser().resolve()
        if path.exists() and not overwrite:
            raise FileExistsError(f"Parquet output already exists: {path}")
        path.parent.mkdir(parents=True, exist_ok=True)

        copy_sql = (
            f"COPY ({self.sql}) TO {_sql_string(str(path))} "
            f"(FORMAT PARQUET, COMPRESSION {_sql_string(compression)})"
        )
        with duckdb.connect(str(self.db_path), read_only=True) as conn:
            conn.execute(copy_sql, list(self.parameters))
        return path


@dataclass(frozen=True)
class GenomeView:
    """Queries for one genome across every sample."""

    db_path: Path
    genome: str

    def genome_stats(self) -> Query:
        """Return genome statistics across samples."""
        return Query(
            self.db_path,
            """
            SELECT sample_id, genome, coverage, breadth, ber
            FROM genome_stats
            WHERE genome = ?
            ORDER BY sample_id
            """,
            (self.genome,),
        )

    def gene_stats(self, gene: str | None = None) -> Query:
        """Return gene statistics across samples."""
        conditions = ["genome = ?"]
        parameters: list[Any] = [self.genome]
        if gene is not None:
            conditions.append("gene = ?")
            parameters.append(_required_text(gene, "gene"))
        return Query(
            self.db_path,
            f"""
            SELECT sample_id, genome, gene, coverage, breadth, ber
            FROM gene_stats
            WHERE {' AND '.join(conditions)}
            ORDER BY sample_id, gene
            """,
            tuple(parameters),
        )

    def profiles(self, gene: str | None = None) -> Query:
        """Return profile positions across samples."""
        conditions = ["genome = ?"]
        parameters: list[Any] = [self.genome]
        if gene is not None:
            conditions.append("gene = ?")
            parameters.append(_required_text(gene, "gene"))
        return Query(
            self.db_path,
            f"""
            SELECT sample_id, chrom, pos, genome, gene, A, C, G, T
            FROM profile_positions
            WHERE {' AND '.join(conditions)}
            """,
            tuple(parameters),
        )

    def sylph_abundance(self) -> Query:
        """Return Sylph rows matching this genome or accession."""
        return Query(
            self.db_path,
            """
            SELECT sample_id, genome, accession, abundance
            FROM sylph_abundance
            WHERE genome = ? OR accession = ?
            ORDER BY sample_id
            """,
            (self.genome, self.genome),
        )


@dataclass(frozen=True)
class SampleView:
    """Queries for one sample across its genomes and genes."""

    db_path: Path
    sample_id: str

    def genome_stats(self, genome: str | None = None) -> Query:
        """Return genome statistics stored for this sample."""
        return self._filtered_query(
            table="genome_stats",
            columns="sample_id, genome, coverage, breadth, ber",
            genome=genome,
            order_by="genome",
        )

    def gene_stats(self, genome: str | None = None, gene: str | None = None) -> Query:
        """Return gene statistics stored for this sample."""
        return self._filtered_query(
            table="gene_stats",
            columns="sample_id, genome, gene, coverage, breadth, ber",
            genome=genome,
            gene=gene,
            order_by="genome, gene",
        )

    def profile(self, genome: str | None = None, gene: str | None = None) -> Query:
        """Return profile positions stored for this sample."""
        return self._filtered_query(
            table="profile_positions",
            columns="sample_id, chrom, pos, genome, gene, A, C, G, T",
            genome=genome,
            gene=gene,
        )

    def sylph_abundance(self, genome: str | None = None) -> Query:
        """Return Sylph abundance rows stored for this sample."""
        conditions = ["sample_id = ?"]
        parameters: list[Any] = [self.sample_id]
        if genome is not None:
            value = _required_text(genome, "genome")
            conditions.append("(genome = ? OR accession = ?)")
            parameters.extend([value, value])
        return Query(
            self.db_path,
            f"""
            SELECT sample_id, genome, accession, abundance
            FROM sylph_abundance
            WHERE {' AND '.join(conditions)}
            ORDER BY genome, accession
            """,
            tuple(parameters),
        )

    def _filtered_query(
        self,
        *,
        table: str,
        columns: str,
        genome: str | None = None,
        gene: str | None = None,
        order_by: str | None = None,
    ) -> Query:
        conditions = ["sample_id = ?"]
        parameters: list[Any] = [self.sample_id]
        if genome is not None:
            conditions.append("genome = ?")
            parameters.append(_required_text(genome, "genome"))
        if gene is not None:
            conditions.append("gene = ?")
            parameters.append(_required_text(gene, "gene"))
        order_clause = f"ORDER BY {order_by}" if order_by else ""
        return Query(
            self.db_path,
            f"""
            SELECT {columns}
            FROM {table}
            WHERE {' AND '.join(conditions)}
            {order_clause}
            """,
            tuple(parameters),
        )


def open_database(db_path: str | Path) -> MetaTrawlDatabase:
    """Open a MetaTrawl database for genome- and sample-centric queries."""
    return MetaTrawlDatabase(db_path)


def _required_text(value: str, name: str) -> str:
    clean = str(value).strip()
    if not clean:
        raise ValueError(f"{name} cannot be empty")
    return clean


def _sql_string(value: str) -> str:
    return "'" + value.replace("'", "''") + "'"

"""Storage-specific implementations of the reference-relative SNP query."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import Iterator

import duckdb
import polars as pl
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq

from metatrawl import allele_mask, provenance


# Full ZipStrain profiles use A/C/G/T bits; the compact format uses A/T/C/G.
_MASK_BITS = dict(zip("ACGT", (1, 4, 8, 2)))
_MASK_BASES = pa.array([next((b for b, bit in _MASK_BITS.items() if bit == n), None)
                        for n in range(16)], type=pa.string())
_MASK_ALTS = pa.array([[b for b, bit in _MASK_BITS.items() if n & bit]
                       for n in range(16)], type=pa.list_(pa.string()))
_MASK_COUNTS = pa.array([n.bit_count() for n in range(16)], type=pa.uint8())
_SNP_SCHEMA = pa.schema([
    ("sample_id", pa.string()),
    ("genome", pa.string()),
    ("chrom", pa.string()),
    ("pos", pa.int64()),
    ("ref_base", pa.string()),
    ("alt_bases", pa.list_(pa.string())),
    ("allele_count", pa.uint8()),
    ("ref_present", pa.bool_()),
])


def full_profile_sql(where: str) -> str:
    """Filter and annotate full profiles inside DuckDB before materialization."""
    bits = tuple(zip("ACGT", (1, 2, 4, 8)))
    alt_bases = ", ".join(
        f"CASE WHEN (alt_mask & {bit}) != 0 THEN ['{base}'] ELSE [] END"
        for base, bit in bits
    )
    allele_mask_sql = " + ".join(
        f"CASE WHEN {base} > 0 THEN {bit} ELSE 0 END" for base, bit in bits
    )
    ref_cases = " ".join(f"WHEN {bit} THEN '{base}'" for base, bit in bits)
    ref_counts = " ".join(f"WHEN {bit} THEN {base}" for base, bit in bits)
    con_cases = " ".join(
        f"WHEN {base} = con_count THEN '{base}'" for base, _ in bits[:-1]
    )
    return f"""
        WITH profile AS (
            SELECT sample_id, genome, chrom, pos, A, C, G, T, ref_base_bitmask,
                   A::BIGINT + C::BIGINT + G::BIGINT + T::BIGINT AS position_coverage,
                   {allele_mask_sql} AS observed_mask
            FROM profile_positions
            WHERE {where}
              AND ref_base_bitmask IN (1, 2, 4, 8)
        ), variants AS (
            SELECT *, observed_mask & ~ref_base_bitmask::INTEGER AS alt_mask,
                   (observed_mask & ref_base_bitmask) != 0 AS ref_present,
                   bit_count(observed_mask)::UTINYINT AS allele_count,
                   CASE ref_base_bitmask {ref_cases} END AS ref_base,
                   CASE ref_base_bitmask {ref_counts} END AS ref_count,
                   greatest(A, C, G, T) AS con_count
            FROM profile
            WHERE position_coverage >= ? AND observed_mask != ref_base_bitmask
        ), consensus AS (
            SELECT *, CASE {con_cases} ELSE 'T' END AS con_base
            FROM variants
        )
        SELECT sample_id, genome, chrom, pos, ref_base,
               list_concat({alt_bases}) AS alt_bases,
               allele_count, ref_present, A, C, G, T, position_coverage, con_base,
               ref_count / position_coverage AS ref_freq,
               con_count / position_coverage AS con_freq,
               (position_coverage - ref_count) / position_coverage AS alt_freq,
               CASE
                   WHEN allele_count = 1 THEN 'SNS'
                   WHEN con_base = ref_base THEN 'SNV'
                   WHEN ref_present THEN 'con_SNV'
                   ELSE 'pop_SNV'
               END AS class
        FROM consensus
    """


@dataclass(frozen=True)
class AlleleMaskSnpQuery:
    """Reference SNP query over compressed blocks, without read-count fields."""

    db_path: Path
    sql: str
    parameters: tuple[str, ...]

    def _batches(self) -> Iterator[pa.RecordBatch]:
        references = allele_mask.ReferenceSegmentCache()
        with duckdb.connect(str(self.db_path), read_only=True) as conn, conn.cursor() as ref_conn:
            provenance.check_schema(conn)
            blocks = conn.execute(self.sql, self.parameters)
            while rows := blocks.fetchmany(16):
                for sample_id, genome, chrom, presence, deviation, covered, checksum in rows:
                    reference = references.segment(ref_conn, genome=genome, chrom=chrom)
                    pos, ref, observed = allele_mask.decode_profile_snps(
                        reference, presence, deviation, covered, checksum, sample_id=sample_id
                    )
                    if not pos.size:
                        continue
                    yield pa.RecordBatch.from_arrays([
                        pa.repeat(sample_id, pos.size),
                        pa.repeat(genome, pos.size),
                        pa.repeat(chrom, pos.size),
                        pa.array(pos, type=pa.int64()),
                        pc.take(_MASK_BASES, pa.array(ref)),
                        pc.take(_MASK_ALTS, pa.array(observed & ~ref)),
                        pc.take(_MASK_COUNTS, pa.array(observed)),
                        pa.array((observed & ref) != 0),
                    ], schema=_SNP_SCHEMA)

    def collect(self) -> pl.DataFrame:
        """Collect variant rows as a Polars DataFrame."""
        return pl.from_arrow(pa.Table.from_batches(self._batches(), schema=_SNP_SCHEMA))

    def lazy(self) -> pl.LazyFrame:
        """Collect variant rows, then return a lazy Polars frame."""
        return self.collect().lazy()

    def sink_parquet(
        self,
        output_file: str | Path,
        *,
        compression: str = "zstd",
        overwrite: bool = False,
    ) -> Path:
        """Write decoded SNPs block by block without collecting them in Python."""
        path = Path(output_file).expanduser().resolve()
        if path.exists() and not overwrite:
            raise FileExistsError(f"Parquet output already exists: {path}")
        path.parent.mkdir(parents=True, exist_ok=True)
        with NamedTemporaryFile(dir=path.parent, prefix=f".{path.name}.", delete=False) as handle:
            staging = Path(handle.name)
        try:
            with pq.ParquetWriter(staging, _SNP_SCHEMA, compression=compression) as writer:
                for batch in self._batches():
                    writer.write_batch(batch)
            if overwrite:
                staging.replace(path)
            else:
                path.hardlink_to(staging)
        finally:
            staging.unlink(missing_ok=True)
        return path

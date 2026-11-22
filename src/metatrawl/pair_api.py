"""On-demand ZipStrain analyses for a single genome and two stored samples."""

from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
import polars as pl
import pyarrow as pa
import pyarrow.parquet as pq

from metatrawl import allele_mask, db


def _manifest(conn, sample: str) -> dict:
    """Recorded provenance for one sample; {} when unknown, including MetaTrawl 0.x files."""
    if not db._table_exists(conn, "sample_provenance"):
        return {}
    row = conn.execute("SELECT manifest_json FROM sample_provenance WHERE sample_id=?", [sample]).fetchone()
    return json.loads(row[0]) if row else {}


@dataclass(frozen=True)
class PairView:
    """A scoped pair; data is fetched only when an analysis is executed."""

    db_path: Path
    sample_1: str
    sample_2: str
    genome: str
    min_cov: int | None = None
    ani_method: str = "popani"
    allow_unknown_provenance: bool = False

    def identical_stretches(self, *, min_length: int = 1) -> PairAnalysisQuery:
        """Return one-based inclusive stretches, with ZipStrain significance estimates."""
        if min_length < 1:
            raise ValueError("min_length must be positive")
        return PairAnalysisQuery(self, "identical_stretches", {"min_length": min_length})

    def mosaic_regions(self, *, factor: float = 20.0, switch_probability: float = 1e-3) -> PairAnalysisQuery:
        """Return ZipStrain's low/background/high-divergence regions on demand."""
        return PairAnalysisQuery(self, "mosaic_regions", {"factor": factor, "switch_probability": switch_probability})

    def compare(self, calculate: str = "genome_ani+max_identical_stretch", *,
                allele_integration: str = "consensus") -> PairAnalysisQuery:
        """Compare genome metrics, or gene selection from retained observed codons."""
        from zipstrain import options
        metrics = options.parse_calculate(calculate)
        if "gene_selection" in metrics:
            if set(metrics) != {"gene_selection"}:
                raise ValueError("Request gene_selection separately from genome metrics in this API")
            return self.gene_selection(allele_integration=allele_integration)
        if not set(metrics) <= set(options.GENOME_METRICS):
            raise ValueError("This pair API supports genome metrics and gene_selection")
        return PairAnalysisQuery(self, "compare", {"calculate": calculate})

    def gene_selection(self, *, allele_integration: str = "consensus", min_cov: int | None = None,
                       min_freq: float = 0.0) -> PairAnalysisQuery:
        """Calculate pairwise dN/dS and pN/pS from retained read-phased codons."""
        if allele_integration not in ("consensus", "weighted"):
            raise ValueError("allele_integration must be consensus or weighted")
        if min_cov is None:
            min_cov = 5 if self.min_cov is None else self.min_cov
        return PairAnalysisQuery(self, "gene_selection", dict(
            allele_integration=allele_integration, min_cov=min_cov, min_freq=min_freq))

    def _validate(self, conn) -> None:
        known = []
        generations = []
        for sample in (self.sample_1, self.sample_2):
            if conn.execute("SELECT 1 FROM samples WHERE sample_id=? AND status='complete'", [sample]).fetchone() is None:
                raise ValueError(f"Unknown or incomplete sample: {sample}")
            manifest = _manifest(conn, sample)
            reference = manifest.get("references", {}).get(self.genome, {})
            contract = manifest.get("contract")
            if (not reference or not contract) and not self.allow_unknown_provenance:
                raise ValueError(f"Unverified provenance for {sample}/{self.genome}; pass allow_unknown_provenance=True to acknowledge it")
            known.append((contract, reference))
            generation = (conn.execute("SELECT profile_schema FROM sample_profile_metadata WHERE sample_id=?", [sample]).fetchone()
                          if db._table_exists(conn, "sample_profile_metadata") else None)
            generations.append(generation[0] if generation else "legacy-unknown")
        if len(set(generations)) != 1:
            raise ValueError("Cannot compare legacy and V2 profiles")
        for index, label in ((0, "profiling contract"), (1, "genome reference/annotation")):
            if all(k[index] for k in known) and known[0][index] != known[1][index]:
                raise ValueError(f"Pair uses incompatible {label}")

    def _execute(self, action):
        from zipstrain.compare import PairComparison

        with TemporaryDirectory(prefix="metatrawl-pair-") as tmp, db.connect_read_only(self.db_path) as conn:
            self._validate(conn)
            storage = db.profile_storage_config(conn)
            min_cov = 5 if self.min_cov is None else self.min_cov
            if storage.mode == "allele-mask":
                if self.ani_method != "popani":
                    raise ValueError("Allele-mask pair analyses support popani only")
                if self.min_cov is not None and self.min_cov != storage.min_cov:
                    raise ValueError(f"This allele-mask database uses fixed min_cov={storage.min_cov}")
                # Values below are allele-presence indicators, not invented read counts.
                # Stored blocks have already passed the fixed import depth threshold.
                min_cov = 1
            profiles = []
            for sample, name in ((self.sample_1, "a"), (self.sample_2, "b")):
                path = Path(tmp) / f"{name}.parquet"
                if storage.mode == "full":
                    quoted = "'" + str(path).replace("'", "''") + "'"
                    conn.execute(f"""
                        COPY (SELECT chrom, genome, pos, A, C, G, T FROM profile_positions
                              WHERE sample_id=? AND genome=? ORDER BY chrom, pos)
                        TO {quoted} (FORMAT PARQUET, COMPRESSION ZSTD)
                    """, [sample, self.genome])
                else:
                    _export_presence(conn, path, sample, self.genome)
                profiles.append(pl.scan_parquet(path))
            pair = PairComparison(*profiles, min_cov=min_cov, ani_method=self.ani_method,
                                  sample_1=self.sample_1, sample_2=self.sample_2)
            return action(pair)


@dataclass(frozen=True)
class PairAnalysisQuery:
    """A pair analysis with the same collect/lazy/sink interface as table queries."""

    pair: PairView
    method: str
    options: dict

    def _execute(self, action):
        if self.method == "gene_selection":
            return _selection(self.pair, self.options, action)
        return self.pair._execute(lambda pair: action(getattr(pair, self.method)(**self.options)))

    def collect(self) -> pl.DataFrame:
        """Execute the scoped analysis and collect its result."""
        return self._execute(lambda frame: frame.collect(engine="streaming"))

    def lazy(self) -> pl.LazyFrame:
        """Execute the analysis, then expose its result as a lazy frame."""
        return self.collect().lazy()

    def sink_parquet(self, output_file, *, compression="zstd", overwrite=False) -> Path:
        """Run the analysis against temporary scoped exports and publish a Parquet."""
        target = Path(output_file).expanduser().resolve()
        if target.exists() and not overwrite:
            raise FileExistsError(f"Parquet output already exists: {target}")
        target.parent.mkdir(parents=True, exist_ok=True)
        with TemporaryDirectory(dir=target.parent) as tmp:
            staging = Path(tmp) / "result.parquet"
            self._execute(lambda frame: frame.sink_parquet(
                staging, compression=compression, engine="streaming"))
            if overwrite:
                staging.replace(target)
            else:
                target.hardlink_to(staging)
        return target


def _export_presence(conn, path: Path, sample: str, genome: str) -> None:
    schema = pa.schema([("chrom", pa.string()), ("genome", pa.string()), ("pos", pa.int64()),
                        *[(base, pa.uint8()) for base in "ACGT"]])
    with pq.ParquetWriter(path, schema, compression="zstd") as writer:
        for block in allele_mask.iter_decoded_sparse_profile_blocks(conn, sample_id=sample, genome=genome):
            values = block.masks
            arrays = [pa.repeat(block.chrom, values.size), pa.repeat(genome, values.size),
                      pa.array(block.positions + block.start_pos)]
            arrays.extend(pa.array(((values & bit) != 0).astype(np.uint8)) for bit in (1, 4, 8, 2))
            writer.write_batch(pa.RecordBatch.from_arrays(arrays, schema=schema))


def _export_query(conn, path: Path, sql: str, parameters: list, metadata: dict) -> None:
    reader = conn.execute(sql, parameters).fetch_record_batch(65_536)
    schema = reader.schema.with_metadata({k.encode(): str(v).encode() for k, v in metadata.items()})
    with pq.ParquetWriter(path, schema, compression="zstd") as writer:
        for batch in reader:
            writer.write_batch(batch)


def _selection(pair: PairView, parameters: dict, action):
    """Export only this genome's authenticated codons; never read profile positions."""
    from zipstrain import codons, gene_info as gi, utils

    with TemporaryDirectory(prefix="metatrawl-selection-") as tmp, db.connect_read_only(pair.db_path) as conn:
        pair._validate(conn)
        if not db._table_exists(conn, "sample_profile_metadata"):
            raise ValueError("Gene selection requires V2 profiles with collect_codons=true")
        evidence = {}
        reference = None
        for sample in (pair.sample_1, pair.sample_2):
            row = conn.execute("SELECT codon_metadata_json FROM sample_profile_metadata WHERE sample_id=?",
                               [sample]).fetchone()
            if row is None or row[0] is None:
                raise ValueError(f"Sample {sample} has no retained codon evidence; profile with collect_codons=true")
            evidence[sample] = json.loads(row[0])
            manifest = _manifest(conn, sample)
            reference = manifest.get("references", {}).get(pair.genome)
            if reference is None:
                raise ValueError("Gene selection requires per-genome reference provenance from the imported bundle")
        if evidence[pair.sample_1].get("zipstrain_codon_filters") != evidence[pair.sample_2].get("zipstrain_codon_filters"):
            raise ValueError("Samples use incompatible codon filters")
        columns = ", ".join(f'"{name}"' for name in gi.GENE_SCHEMA)
        genes = conn.execute(f"SELECT {columns} FROM gene_catalog WHERE genome=? ORDER BY gene_id",
                             [pair.genome]).pl().cast(gi.GENE_SCHEMA)
        if genes.is_empty():
            raise ValueError(f"No gene catalog for genome {pair.genome}")
        # Authenticate a genome-scoped projection of the imported reference and
        # catalog. IDs stay unchanged; parent sample receipts remain in DuckDB.
        identity = gi.content_digest(sorted(
            (chrom, pair.genome, length, digest) for chrom, length, digest in reference["scaffolds"]))
        meta = {gi.SCHEMA_KEY: "2", gi.REFERENCE_ID: identity, gi.ANNOTATION_ID: gi.content_digest(genes.rows())}
        root = Path(tmp)
        _export_query(conn, root/"reference_codons.parquet", """
            SELECT r.gene_id, r.codon_index, r.reference_codon
            FROM reference_codons r JOIN gene_catalog g USING(gene_id)
            WHERE g.genome=? ORDER BY r.gene_id, r.codon_index
        """, [pair.genome], meta)
        meta["zipstrain_reference_codons_sha256"] = utils.sha256_file(root/"reference_codons.parquet")
        catalog = root/"gene_info.parquet"
        genes.write_parquet(catalog, compression="zstd", metadata=meta)
        paths = []
        for index, sample in enumerate((pair.sample_1, pair.sample_2)):
            path = root/f"sample{index}.parquet"
            _export_query(conn, path, """
                SELECT c.gene_id, c.codon_index, c.codon, c.count
                FROM codon_counts c JOIN gene_catalog g USING(gene_id)
                WHERE c.sample_id=? AND g.genome=? ORDER BY c.gene_id, c.codon_index, c.codon
            """, [sample, pair.genome], {**evidence[sample], **meta})
            paths.append(path)
        rates = codons.pairwise_selection(*paths, catalog, **parameters)
        result = rates.join(genes.lazy().select("gene_id", "gene", "genome", "chrom", "start", "end"),
                            on="gene_id", how="left").with_columns(
                                sample_1=pl.lit(pair.sample_1), sample_2=pl.lit(pair.sample_2))
        return action(result)

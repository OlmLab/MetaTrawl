"""Import authenticated ZipStrain V2 catalogs and read-phased codon evidence."""

from __future__ import annotations

import json
import shutil
from pathlib import Path
from tempfile import TemporaryDirectory

import polars as pl

from zipstrain import codons, gene_info as gi, utils


def metadata(path: Path) -> dict[str, str]:
    """Read scientific footer fields without retaining the Arrow schema blob."""
    return {k: v for k, v in pl.read_parquet_metadata(path).items() if k.startswith("zipstrain_")}


def import_evidence(conn, bundle) -> None:
    """Store catalog/evidence inside the caller's sample import transaction.

    Only the position footer is read here, never the database position table.
    Shared reference codons are inserted only for previously unseen gene IDs.
    """
    profile_meta = metadata(bundle.profile_file)
    generation = profile_meta.get(gi.SCHEMA_KEY, "legacy-unknown")
    if generation not in ("legacy-unknown", "2"):
        raise ValueError(f"sample={bundle.run_id} unsupported ZipStrain profile schema: {generation}")
    counts_meta = None
    if generation == "2":
        if bundle.gene_info_file is None or bundle.reference_codons_file is None:
            raise ValueError(f"sample={bundle.run_id} V2 import requires gene_info_file and reference_codons_file")
        with TemporaryDirectory(prefix="metatrawl-v2-assets-") as tmp:
            catalog = Path(tmp) / "gene_info.parquet"
            reference = Path(tmp) / "reference_codons.parquet"
            # ZipStrain resolves the catalog path before finding its sibling
            # reference codons, so a catalog symlink would lose this layout.
            shutil.copyfile(bundle.gene_info_file, catalog)
            reference.symlink_to(Path(bundle.reference_codons_file).resolve())
            gi.validate_profile_assets([bundle.profile_file, reference], catalog)
            expected = metadata(catalog)
            if bundle.gene_stats_file is not None:
                gi.validate_profile_assets([bundle.gene_stats_file], catalog)
                stats_meta = metadata(bundle.gene_stats_file)
                if (stats_meta.get("zipstrain_selection_method", "not_calculated") != "not_calculated"
                        and bundle.codon_counts_file is None):
                    raise ValueError(f"sample={bundle.run_id} selection statistics require the retained codon-count table")
            if utils.sha256_file(reference) != expected.get("zipstrain_reference_codons_sha256"):
                raise ValueError("Reference codon contents do not match the gene catalog")
            columns = list(gi.GENE_SCHEMA)
            quoted = ", ".join(f'"{name}"' for name in columns)
            # The catalog was authenticated above. A native Parquet scan avoids
            # signed-integer coercion in Arrow filter pushdown for UInt64 IDs.
            conn.execute(f"""
                CREATE OR REPLACE TEMP TABLE _v2_genes AS
                SELECT {quoted} FROM read_parquet(?)
            """, [str(catalog)])
            if bundle.gene_stats_file is not None:
                if conn.execute("""
                    SELECT EXISTS(
                        SELECT 1 FROM read_parquet(?) s LEFT JOIN _v2_genes g USING(gene_id)
                        WHERE g.gene_id IS NULL OR s.genome IS DISTINCT FROM g.genome
                            OR s.chrom IS DISTINCT FROM g.chrom OR s.start IS DISTINCT FROM g.start
                            OR s.end IS DISTINCT FROM g.end OR s.length IS DISTINCT FROM g.length
                    )
                """, [str(bundle.gene_stats_file)]).fetchone()[0]:
                    raise ValueError(f"sample={bundle.run_id} gene statistics do not match the gene catalog")
            tests = " OR ".join(f'c."{name}" IS DISTINCT FROM n."{name}"' for name in columns if name != "gene_id")
            if conn.execute(f"""
                SELECT EXISTS(SELECT 1 FROM gene_catalog c JOIN _v2_genes n USING(gene_id)
                WHERE {tests})
            """).fetchone()[0]:
                raise ValueError("Gene ID collision or incompatible gene catalog")
            conn.execute("""
                CREATE OR REPLACE TEMP TABLE _v2_new_gene_ids AS
                SELECT gene_id FROM _v2_genes n
                WHERE NOT EXISTS (SELECT 1 FROM gene_catalog c WHERE c.gene_id=n.gene_id)
            """)
            conn.execute("""
                INSERT INTO reference_codons
                SELECT r.gene_id, r.codon_index, r.reference_codon FROM read_parquet(?) r
                JOIN _v2_new_gene_ids n USING (gene_id)
            """, [str(reference)])
            conn.execute(f"""
                INSERT INTO gene_catalog ({quoted})
                SELECT {quoted} FROM _v2_genes
                WHERE gene_id IN (SELECT gene_id FROM _v2_new_gene_ids)
            """)
            conn.execute("DROP TABLE IF EXISTS _v2_new_gene_ids")
            conn.execute("DROP TABLE _v2_genes")
            conn.execute("INSERT OR IGNORE INTO annotation_catalogs VALUES (?, ?, ?)",
                         [expected[gi.ANNOTATION_ID], expected[gi.REFERENCE_ID], json.dumps(expected, sort_keys=True)])
            if bundle.codon_counts_file is not None:
                counts_meta = metadata(bundle.codon_counts_file)
                codons.read_codon_counts(bundle.codon_counts_file, catalog)
                if counts_meta.get("zipstrain_source_profile_sha256") != utils.sha256_file(bundle.profile_file):
                    raise ValueError(f"sample={bundle.run_id} codon evidence does not match its position profile")
                conn.execute("""
                    INSERT INTO codon_counts
                    SELECT ?, gene_id::UBIGINT, codon_index::UINTEGER, codon::VARCHAR, count::UBIGINT
                    FROM read_parquet(?)
                """, [bundle.run_id, str(bundle.codon_counts_file)])
    elif bundle.codon_counts_file is not None:
        raise ValueError("Codon evidence requires a V2 position profile")
    conn.execute("INSERT OR REPLACE INTO sample_profile_metadata VALUES (?, ?, ?, ?, ?)",
                 [bundle.run_id, generation, json.dumps(profile_meta, sort_keys=True),
                  json.dumps(counts_meta, sort_keys=True) if counts_meta is not None else None,
                  counts_meta is not None])

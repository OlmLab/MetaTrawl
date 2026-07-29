from __future__ import annotations

from pathlib import Path

import duckdb
import h5py
import numpy as np
import polars as pl
import pytest

from metatrawl import allele_mask
from metatrawl import db
from metatrawl import matrix_hdf5
from metatrawl import migration
from metatrawl.api import ProfileCountsUnavailableError, open_database


GENOME = "GCF_000001.1"
SAMPLES = ("S1", "S2", "S3")


def _project_files(tmp_path: Path) -> tuple[Path, Path, Path, Path]:
    cache_dir = tmp_path / "cache"
    (cache_dir / "genomes").mkdir(parents=True)
    (cache_dir / "genomes" / f"{GENOME}.fna").write_text(
        ">contigA\nATCGNATCGT\n"
        ">contigB\nGGCAT\n"
    )
    bed_file = tmp_path / "genome.bed"
    bed_file.write_text("contigA\t0\t10\ncontigB\t0\t5\n")
    stb_file = tmp_path / "genome.stb"
    stb_file.write_text(f"contigA\t{GENOME}\ncontigB\t{GENOME}\n")
    gene_ranges = tmp_path / "genes.tsv"
    gene_ranges.write_text(
        "gene_a\tcontigA\t1\t6\n"
        "gene_b\tcontigB\t1\t5\n"
    )
    return cache_dir, bed_file, stb_file, gene_ranges


def _profile_rows(sample: str) -> list[tuple[str, int, int, int, int, int]]:
    rows = {
        "S1": [
            ("contigA", 1, 5, 0, 0, 0),   # A
            ("contigA", 2, 0, 0, 0, 5),   # T
            ("contigA", 3, 0, 3, 3, 0),   # C+G
            ("contigA", 4, 0, 0, 5, 0),   # G
            ("contigA", 5, 5, 0, 0, 0),   # ambiguous reference
            ("contigA", 9, 0, 0, 0, 4),   # below threshold
            ("contigB", 1, 0, 0, 5, 0),
            ("contigB", 2, 0, 0, 4, 0),   # below threshold
            ("contigB", 3, 0, 5, 0, 0),
        ],
        "S2": [
            ("contigA", 1, 5, 0, 0, 0),
            ("contigA", 2, 0, 0, 0, 5),
            ("contigA", 3, 0, 5, 0, 0),
            ("contigA", 4, 5, 0, 0, 0),   # no shared allele with S1
            ("contigA", 5, 0, 0, 0, 5),   # ambiguous reference
            ("contigB", 1, 0, 0, 5, 0),
            ("contigB", 3, 0, 0, 0, 5),   # no shared allele with S1
        ],
        "S3": [
            ("contigA", 1, 4, 0, 0, 0),
            ("contigB", 1, 0, 0, 4, 0),
        ],
    }
    return rows[sample]


def _write_bundle(tmp_path: Path, sample: str) -> db.ProfileBundle:
    rows = _profile_rows(sample)
    profile_file = tmp_path / f"{sample}.profile.parquet"
    pl.DataFrame(
        {
            "chrom": [row[0] for row in rows],
            "genome": [GENOME] * len(rows),
            "gene": ["gene_a" if row[0] == "contigA" else "gene_b" for row in rows],
            "pos": [row[1] for row in rows],
            "A": [row[2] for row in rows],
            "C": [row[3] for row in rows],
            "G": [row[4] for row in rows],
            "T": [row[5] for row in rows],
        }
    ).write_parquet(profile_file)
    genome_stats = tmp_path / f"{sample}.genome_stats.parquet"
    pl.DataFrame(
        {
            "genome": [GENOME],
            "coverage": [5.0],
            "breadth": [0.8],
            "ber": [0.9],
            "ref_ani": [0.999],
        }
    ).write_parquet(genome_stats)
    gene_stats = tmp_path / f"{sample}.gene_stats.parquet"
    pl.DataFrame(
        {
            "genome": [GENOME, GENOME],
            "gene": ["gene_a", "gene_b"],
            "coverage": [5.0, 5.0],
            "breadth": [0.8, 0.8],
            "ber": [0.9, 0.9],
            "ref_ani": [0.999, 0.999],
            "length": [6, 5],
        }
    ).write_parquet(gene_stats)
    sylph = tmp_path / f"{sample}.sylph.csv"
    pl.DataFrame(
        {
            "genome": [GENOME],
            "accession": [GENOME],
            "abundance": [0.2],
        }
    ).write_csv(sylph)
    return db.ProfileBundle(
        run_id=sample,
        profile_file=profile_file,
        genome_stats_file=genome_stats,
        gene_stats_file=gene_stats,
        sylph_abundance_file=sylph,
    )


def _make_database(
    tmp_path: Path,
    *,
    name: str,
    mode: str,
    cache_dir: Path,
    samples: tuple[str, ...] = SAMPLES,
) -> Path:
    db_file = tmp_path / f"{name}.duckdb"
    with db.connect(db_file) as conn:
        db.configure_profile_storage(
            conn,
            mode=mode,
            min_cov=5 if mode == "allele-mask" else None,
        )
        for sample in samples:
            db.add_runs(conn, [sample])
            db.import_profile_bundle(
                conn,
                _write_bundle(tmp_path, sample),
                cache_dir=cache_dir,
            )
    return db_file


def _build_matrix(
    db_file: Path,
    output_file: Path,
    *,
    bed_file: Path,
    stb_file: Path,
    gene_ranges: Path,
    sparse: bool,
) -> None:
    with db.connect(db_file) as conn:
        matrix_hdf5.build_matrix_hdf5_from_duckdb(
            conn,
            sample_ids=list(SAMPLES),
            output_file=output_file,
            genome=GENOME,
            bed_file=bed_file,
            stb_file=stb_file,
            gene_range_table=gene_ranges,
            min_cov=5,
            sparse=sparse,
        )


def _hdf_matrix_payload(path: Path) -> tuple:
    with h5py.File(path, "r") as handle:
        matrix = handle["matrices"]["0"]
        if isinstance(matrix, h5py.Group):
            return (
                matrix["indptr"][...].tolist(),
                matrix["indices"][...].tolist(),
                matrix["values"][...].tolist(),
            )
        return (matrix[...].tolist(),)


def _compare_results(path: Path) -> tuple[list[tuple], list[tuple]]:
    with duckdb.connect(str(path), read_only=True) as conn:
        genome_rows = conn.execute(
            """
            SELECT sample_1, sample_2, genome, total_positions,
                   share_allele_pos, genome_ani, max_consecutive_length
            FROM matrix_compare_results
            ORDER BY sample_idx_1, sample_idx_2, genome_idx
            """
        ).fetchall()
        gene_rows = conn.execute(
            """
            SELECT sample_1, sample_2, genome, gene, gene_pop_ani
            FROM matrix_compare_gene_results
            ORDER BY sample_idx_1, sample_idx_2, genome_idx, gene
            """
        ).fetchall()
    return genome_rows, gene_rows


def test_nibble_and_presence_codecs_round_trip_odd_lengths() -> None:
    masks = np.asarray([0, 1, 2, 4, 8, 15, 3], dtype=np.uint8)
    presence = masks > 0
    assert np.array_equal(
        allele_mask.unpack_nibbles(allele_mask.pack_nibbles(masks), len(masks)),
        masks,
    )
    assert np.array_equal(
        allele_mask.unpack_presence(
            allele_mask.pack_presence(presence),
            len(presence),
        ),
        presence,
    )


def test_allele_mask_import_keeps_stats_but_not_profile_rows(tmp_path: Path) -> None:
    cache_dir, *_ = _project_files(tmp_path)
    compact = _make_database(
        tmp_path,
        name="compact",
        mode="allele-mask",
        cache_dir=cache_dir,
    )
    with db.connect(compact) as conn:
        storage = db.profile_storage_config(conn)
        assert storage.mode == "allele-mask"
        assert storage.min_cov == 5
        assert conn.execute("SELECT count(*) FROM profile_positions").fetchone() == (0,)
        assert conn.execute("SELECT count(*) FROM genome_stats").fetchone() == (3,)
        assert conn.execute("SELECT count(*) FROM gene_stats").fetchone() == (6,)
        assert conn.execute("SELECT count(*) FROM sylph_abundance").fetchone() == (3,)
        assert conn.execute(
            "SELECT count(*) FROM allele_mask_reference_segments"
        ).fetchone() == (2,)
        # S3 has only sub-threshold rows and is represented by an implicit zero row.
        assert conn.execute(
            "SELECT count(*) FROM allele_mask_profile_blocks WHERE sample_id = 'S3'"
        ).fetchone() == (0,)


@pytest.mark.parametrize("sparse", [False, True])
def test_full_and_allele_mask_hdf5_payloads_are_identical(
    tmp_path: Path,
    sparse: bool,
) -> None:
    cache_dir, bed_file, stb_file, gene_ranges = _project_files(tmp_path)
    full = _make_database(tmp_path, name="full", mode="full", cache_dir=cache_dir)
    compact = _make_database(
        tmp_path,
        name="compact",
        mode="allele-mask",
        cache_dir=cache_dir,
    )
    full_h5 = tmp_path / f"full-{sparse}.h5"
    compact_h5 = tmp_path / f"compact-{sparse}.h5"
    _build_matrix(
        full,
        full_h5,
        bed_file=bed_file,
        stb_file=stb_file,
        gene_ranges=gene_ranges,
        sparse=sparse,
    )
    _build_matrix(
        compact,
        compact_h5,
        bed_file=bed_file,
        stb_file=stb_file,
        gene_ranges=gene_ranges,
        sparse=sparse,
    )
    assert _hdf_matrix_payload(full_h5) == _hdf_matrix_payload(compact_h5)
    with h5py.File(compact_h5, "r") as handle:
        assert handle["metadata"].attrs["storage_mode"] == "bitmask"
        assert handle["metadata"].attrs["coverage_filter_min_cov"] == "5"
        assert handle["metadata"].attrs["profile_format"] == "metatrawl_allele_mask_v1"
        assert handle["samples"]["sample_name"].asstr()[...].tolist() == list(SAMPLES)


def test_full_and_allele_mask_popani_ibs_and_gene_results_are_identical(
    tmp_path: Path,
) -> None:
    from zipstrain import matrix_pairs

    cache_dir, bed_file, stb_file, gene_ranges = _project_files(tmp_path)
    full = _make_database(tmp_path, name="full", mode="full", cache_dir=cache_dir)
    compact = _make_database(
        tmp_path,
        name="compact",
        mode="allele-mask",
        cache_dir=cache_dir,
    )
    full_h5 = tmp_path / "full.h5"
    compact_h5 = tmp_path / "compact.h5"
    _build_matrix(
        full,
        full_h5,
        bed_file=bed_file,
        stb_file=stb_file,
        gene_ranges=gene_ranges,
        sparse=True,
    )
    _build_matrix(
        compact,
        compact_h5,
        bed_file=bed_file,
        stb_file=stb_file,
        gene_ranges=gene_ranges,
        sparse=True,
    )
    full_compare = tmp_path / "full-compare.duckdb"
    compact_compare = tmp_path / "compact-compare.duckdb"
    for matrix_file, compare_file in (
        (full_h5, full_compare),
        (compact_h5, compact_compare),
    ):
        matrix_pairs.matrix_compare(
            matrix_db_file=matrix_file,
            output_file=compare_file,
            calculate="all",
            ani_method="popani",
            backend="numpy",
            min_cov=5,
        )
    assert _compare_results(full_compare) == _compare_results(compact_compare)


def test_allele_mask_matrix_rejects_counts_and_threshold_changes(
    tmp_path: Path,
) -> None:
    cache_dir, bed_file, stb_file, _ = _project_files(tmp_path)
    compact = _make_database(
        tmp_path,
        name="compact",
        mode="allele-mask",
        cache_dir=cache_dir,
        samples=("S1",),
    )
    with db.connect(compact) as conn:
        with pytest.raises(ValueError, match="only build bitmask"):
            matrix_hdf5.build_matrix_hdf5_from_duckdb(
                conn,
                sample_ids=["S1"],
                output_file=tmp_path / "counts.h5",
                genome=GENOME,
                bed_file=bed_file,
                stb_file=stb_file,
                storage_mode="counts",
                min_cov=5,
            )
        with pytest.raises(ValueError, match="stored with min_cov=5"):
            matrix_hdf5.build_matrix_hdf5_from_duckdb(
                conn,
                sample_ids=["S1"],
                output_file=tmp_path / "wrong-threshold.h5",
                genome=GENOME,
                bed_file=bed_file,
                stb_file=stb_file,
                min_cov=4,
            )
    assert not (tmp_path / "counts.h5").exists()
    assert not (tmp_path / "wrong-threshold.h5").exists()


def test_allele_mask_append_matches_fresh_full_matrix(tmp_path: Path) -> None:
    cache_dir, bed_file, stb_file, gene_ranges = _project_files(tmp_path)
    compact = _make_database(
        tmp_path,
        name="compact",
        mode="allele-mask",
        cache_dir=cache_dir,
        samples=("S1",),
    )
    appended_h5 = tmp_path / "appended.h5"
    with db.connect(compact) as conn:
        matrix_hdf5.build_matrix_hdf5_from_duckdb(
            conn,
            sample_ids=["S1"],
            output_file=appended_h5,
            genome=GENOME,
            bed_file=bed_file,
            stb_file=stb_file,
            gene_range_table=gene_ranges,
            min_cov=5,
            sparse=True,
        )
        for sample in ("S2", "S3"):
            db.add_runs(conn, [sample])
            db.import_profile_bundle(
                conn,
                _write_bundle(tmp_path, sample),
                cache_dir=cache_dir,
            )
        matrix_hdf5.append_matrix_hdf5_from_duckdb(
            conn,
            sample_ids=["S2", "S3"],
            matrix_hdf5_file=appended_h5,
        )

    full = _make_database(tmp_path, name="full", mode="full", cache_dir=cache_dir)
    fresh_h5 = tmp_path / "fresh.h5"
    _build_matrix(
        full,
        fresh_h5,
        bed_file=bed_file,
        stb_file=stb_file,
        gene_ranges=gene_ranges,
        sparse=True,
    )
    assert _hdf_matrix_payload(appended_h5) == _hdf_matrix_payload(fresh_h5)


def test_allele_mask_profile_api_rejects_count_queries(tmp_path: Path) -> None:
    cache_dir, *_ = _project_files(tmp_path)
    compact = _make_database(
        tmp_path,
        name="compact",
        mode="allele-mask",
        cache_dir=cache_dir,
        samples=("S1",),
    )
    database = open_database(compact)
    assert database.sample("S1").genome_stats().collect().height == 1
    with pytest.raises(ProfileCountsUnavailableError, match="does not retain A/C/G/T"):
        database.sample("S1").profile()
    with pytest.raises(ProfileCountsUnavailableError, match="does not retain A/C/G/T"):
        database.genome(GENOME).profiles()


def test_allele_mask_import_rolls_back_without_cached_reference(tmp_path: Path) -> None:
    db_file = tmp_path / "compact.duckdb"
    with db.connect(db_file) as conn:
        db.configure_profile_storage(conn, mode="allele-mask", min_cov=5)
        db.add_runs(conn, ["S1"])
        with pytest.raises(FileNotFoundError, match="cached FASTA"):
            db.import_profile_bundle(conn, _write_bundle(tmp_path, "S1"))
        assert conn.execute("SELECT count(*) FROM samples").fetchone() == (0,)
        assert conn.execute(
            "SELECT count(*) FROM allele_mask_profile_blocks"
        ).fetchone() == (0,)
        assert db.remaining_runs(conn) == ["S1"]


def test_profile_storage_cannot_change_after_import(tmp_path: Path) -> None:
    cache_dir, *_ = _project_files(tmp_path)
    compact = _make_database(
        tmp_path,
        name="compact",
        mode="allele-mask",
        cache_dir=cache_dir,
        samples=("S1",),
    )
    with db.connect(compact) as conn:
        with pytest.raises(ValueError, match="Cannot change profile storage"):
            db.configure_profile_storage(conn, mode="full")
        with pytest.raises(ValueError, match="Cannot change profile storage"):
            db.configure_profile_storage(conn, mode="allele-mask", min_cov=6)


def test_full_database_migration_is_resumable_and_matrix_equivalent(
    tmp_path: Path,
) -> None:
    cache_dir, bed_file, stb_file, gene_ranges = _project_files(tmp_path)
    source = _make_database(tmp_path, name="source", mode="full", cache_dir=cache_dir)
    target = tmp_path / "migrated.duckdb"
    first = migration.migrate_full_database(
        source_db=source,
        output_db=target,
        cache_dir=cache_dir,
        min_cov=5,
    )
    assert first.completed_samples == 3
    assert first.failed_samples == 0
    resumed = migration.migrate_full_database(
        source_db=source,
        output_db=target,
        cache_dir=cache_dir,
        min_cov=5,
    )
    assert resumed.migrated_samples == 0
    assert resumed.completed_samples == 3

    source_h5 = tmp_path / "source.h5"
    target_h5 = tmp_path / "target.h5"
    _build_matrix(
        source,
        source_h5,
        bed_file=bed_file,
        stb_file=stb_file,
        gene_ranges=gene_ranges,
        sparse=True,
    )
    _build_matrix(
        target,
        target_h5,
        bed_file=bed_file,
        stb_file=stb_file,
        gene_ranges=gene_ranges,
        sparse=True,
    )
    assert _hdf_matrix_payload(source_h5) == _hdf_matrix_payload(target_h5)


def test_full_database_migration_recovers_incomplete_initial_copy(
    tmp_path: Path,
) -> None:
    cache_dir, *_ = _project_files(tmp_path)
    source = _make_database(
        tmp_path,
        name="source",
        mode="full",
        cache_dir=cache_dir,
        samples=("S1",),
    )
    target = tmp_path / "interrupted.duckdb"
    # A killed initialization transaction leaves only the empty schema.
    with db.connect(target):
        pass

    result = migration.migrate_full_database(
        source_db=source,
        output_db=target,
        cache_dir=cache_dir,
        min_cov=5,
    )

    assert result.completed_samples == 1
    with db.connect(target) as conn:
        assert conn.execute(
            "SELECT run_id FROM sra_runs ORDER BY run_id"
        ).fetchall() == [("S1",)]
        assert conn.execute(
            "SELECT status FROM allele_mask_migration_state WHERE sample_id = 'S1'"
        ).fetchone() == ("done",)

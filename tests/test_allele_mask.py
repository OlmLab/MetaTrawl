from __future__ import annotations

from pathlib import Path
import errno
import signal

import duckdb
import h5py
import numpy as np
import polars as pl
import pyarrow as pa
import pytest

from metatrawl import allele_mask
from metatrawl import db
from metatrawl import matrix_hdf5
from metatrawl import migration
from metatrawl import workflows
from metatrawl.api import ProfileCountsUnavailableError, open_database
from metatrawl.logging import ThrottledMatrixLogger


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
    positions = np.asarray([0, 3, 6], dtype=np.int64)
    assert np.array_equal(
        allele_mask.unpack_nibbles_at_positions(
            allele_mask.pack_nibbles(masks),
            positions,
        ),
        masks[positions],
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
        if not sparse:
            assert handle["matrices"]["0"].chunks[0] == 1


def test_allele_mask_matrix_reuses_decoded_reference_segments(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    cache_dir, *_ = _project_files(tmp_path)
    compact = _make_database(tmp_path, name="reference-cache", mode="allele-mask", cache_dir=cache_dir)
    original_loader = allele_mask._load_reference_segments
    reference_loads = 0

    def count_reference_loads(*args, **kwargs):
        nonlocal reference_loads
        reference_loads += 1
        return original_loader(*args, **kwargs)

    monkeypatch.setattr(allele_mask, "_load_reference_segments", count_reference_loads)
    reference_cache: dict[str, dict[int, allele_mask.ReferenceSegment]] = {}
    with db.connect(compact) as conn:
        list(
            allele_mask.iter_decoded_profile_blocks(
                conn,
                sample_id="S1",
                genome=GENOME,
                reference_cache=reference_cache,
            )
        )
        list(
            allele_mask.iter_decoded_profile_blocks(
                conn,
                sample_id="S2",
                genome=GENOME,
                reference_cache=reference_cache,
            )
        )
    assert reference_loads == 1


@pytest.mark.parametrize("sparse", [False, True])
def test_matrix_build_resumes_after_a_committed_sample_batch(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    sparse: bool,
) -> None:
    cache_dir, bed_file, stb_file, gene_ranges = _project_files(tmp_path)
    db_file = _make_database(tmp_path, name="resumable", mode="full", cache_dir=cache_dir)
    resumed_h5 = tmp_path / f"resumed-{sparse}.h5"
    checkpoint_h5 = resumed_h5.with_suffix(resumed_h5.suffix + ".tmp")
    original_loader = matrix_hdf5._load_duckdb_sample_batch_genome_matrices

    def interrupt_second_batch(*args, **kwargs):
        if "S2" in kwargs["sample_ids"]:
            raise RuntimeError("simulated cancellation")
        return original_loader(*args, **kwargs)

    monkeypatch.setattr(
        matrix_hdf5,
        "_load_duckdb_sample_batch_genome_matrices",
        interrupt_second_batch,
    )
    with db.connect(db_file) as conn:
        with pytest.raises(RuntimeError, match="simulated cancellation"):
            matrix_hdf5.build_matrix_hdf5_from_duckdb(
                conn,
                sample_ids=list(SAMPLES),
                output_file=resumed_h5,
                genome=GENOME,
                bed_file=bed_file,
                stb_file=stb_file,
                gene_range_table=gene_ranges,
                min_cov=5,
                export_batch_mb=0.000001,
                sparse=sparse,
            )

    assert checkpoint_h5.exists()
    assert not resumed_h5.exists()
    with h5py.File(checkpoint_h5, "r") as handle:
        assert handle.attrs[matrix_hdf5.MATRIX_BUILD_STATE_ATTR] == "building"
        assert handle.attrs[matrix_hdf5.MATRIX_COMMITTED_SAMPLE_COUNT_ATTR] == 1
        assert handle["samples"]["sample_name"].asstr()[...].tolist() == ["S1"]

    monkeypatch.setattr(
        matrix_hdf5,
        "_load_duckdb_sample_batch_genome_matrices",
        original_loader,
    )
    resume_events: list[dict[str, object]] = []
    with db.connect(db_file) as conn:
        matrix_hdf5.build_matrix_hdf5_from_duckdb(
            conn,
            sample_ids=list(SAMPLES),
            output_file=resumed_h5,
            genome=GENOME,
            bed_file=bed_file,
            stb_file=stb_file,
            gene_range_table=gene_ranges,
            min_cov=5,
            export_batch_mb=0.000001,
            sparse=sparse,
            progress_callback=resume_events.append,
        )

    clean_h5 = tmp_path / f"clean-{sparse}.h5"
    _build_matrix(
        db_file,
        clean_h5,
        bed_file=bed_file,
        stb_file=stb_file,
        gene_ranges=gene_ranges,
        sparse=sparse,
    )
    assert not checkpoint_h5.exists()
    assert any(event["phase"] == "resume" for event in resume_events)
    assert _hdf_matrix_payload(resumed_h5) == _hdf_matrix_payload(clean_h5)
    with h5py.File(resumed_h5, "r") as handle:
        assert handle.attrs[matrix_hdf5.MATRIX_BUILD_STATE_ATTR] == "complete"
        assert handle.attrs[matrix_hdf5.MATRIX_COMMITTED_SAMPLE_COUNT_ATTR] == len(SAMPLES)
        assert handle["samples"]["sample_name"].asstr()[...].tolist() == list(SAMPLES)


def test_matrix_write_sigterm_guard_preserves_normal_handler() -> None:
    previous_handler = signal.getsignal(signal.SIGTERM)
    checkpoint = Path("matrix.h5.tmp")
    with pytest.raises(RuntimeError, match="Rerun the same command to resume"):
        with matrix_hdf5._matrix_write_termination_guard(checkpoint):
            active_handler = signal.getsignal(signal.SIGTERM)
            assert callable(active_handler)
            active_handler(signal.SIGTERM, None)
    assert signal.getsignal(signal.SIGTERM) == previous_handler


@pytest.mark.parametrize("sparse", [False, True])
def test_matrix_build_discards_an_uncommitted_partial_batch_on_resume(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    sparse: bool,
) -> None:
    cache_dir, bed_file, stb_file, gene_ranges = _project_files(tmp_path)
    db_file = _make_database(tmp_path, name="partial", mode="full", cache_dir=cache_dir)
    resumed_h5 = tmp_path / f"partial-{sparse}.h5"
    checkpoint_h5 = resumed_h5.with_suffix(resumed_h5.suffix + ".tmp")
    original_loader = matrix_hdf5._load_duckdb_sample_batch_genome_matrices

    def interrupt_batch(*args, **kwargs):
        if "S2" in kwargs["sample_ids"]:
            raise RuntimeError("simulated cancellation")
        return original_loader(*args, **kwargs)

    monkeypatch.setattr(
        matrix_hdf5,
        "_load_duckdb_sample_batch_genome_matrices",
        interrupt_batch,
    )
    with db.connect(db_file) as conn:
        with pytest.raises(RuntimeError, match="simulated cancellation"):
            matrix_hdf5.build_matrix_hdf5_from_duckdb(
                conn,
                sample_ids=list(SAMPLES),
                output_file=resumed_h5,
                genome=GENOME,
                bed_file=bed_file,
                stb_file=stb_file,
                gene_range_table=gene_ranges,
                min_cov=5,
                export_batch_mb=128,
                sparse=sparse,
            )

    with h5py.File(checkpoint_h5, "r") as handle:
        assert handle.attrs[matrix_hdf5.MATRIX_COMMITTED_SAMPLE_COUNT_ATTR] == 0
        assert handle["samples"]["sample_name"].shape == (0,)

    monkeypatch.setattr(
        matrix_hdf5,
        "_load_duckdb_sample_batch_genome_matrices",
        original_loader,
    )
    with db.connect(db_file) as conn:
        matrix_hdf5.build_matrix_hdf5_from_duckdb(
            conn,
            sample_ids=list(SAMPLES),
            output_file=resumed_h5,
            genome=GENOME,
            bed_file=bed_file,
            stb_file=stb_file,
            gene_range_table=gene_ranges,
            min_cov=5,
            export_batch_mb=128,
            sparse=sparse,
        )

    clean_h5 = tmp_path / f"partial-clean-{sparse}.h5"
    _build_matrix(
        db_file,
        clean_h5,
        bed_file=bed_file,
        stb_file=stb_file,
        gene_ranges=gene_ranges,
        sparse=sparse,
    )
    assert _hdf_matrix_payload(resumed_h5) == _hdf_matrix_payload(clean_h5)


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


@pytest.mark.parametrize("sparse", [False, True])
def test_matrix_append_resumes_from_committed_batches(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    sparse: bool,
) -> None:
    cache_dir, bed_file, stb_file, gene_ranges = _project_files(tmp_path)
    compact = _make_database(
        tmp_path,
        name=f"append-resume-{sparse}",
        mode="allele-mask",
        cache_dir=cache_dir,
        samples=("S1",),
    )
    matrix_file = tmp_path / f"append-resume-{sparse}.h5"
    with db.connect(compact) as conn:
        matrix_hdf5.build_matrix_hdf5_from_duckdb(
            conn,
            sample_ids=["S1"],
            output_file=matrix_file,
            genome=GENOME,
            bed_file=bed_file,
            stb_file=stb_file,
            gene_range_table=gene_ranges,
            min_cov=5,
            sparse=sparse,
        )
        for sample in ("S2", "S3"):
            db.add_runs(conn, [sample])
            db.import_profile_bundle(conn, _write_bundle(tmp_path, sample), cache_dir=cache_dir)

        if sparse:
            original_loader = matrix_hdf5._load_allele_mask_sample_genome_sparse

            def interrupt_third_batch(*args, **kwargs):
                if kwargs["sample_id"] == "S3":
                    raise RuntimeError("simulated append cancellation")
                return original_loader(*args, **kwargs)

            loader_name = "_load_allele_mask_sample_genome_sparse"
        else:
            original_loader = matrix_hdf5._load_duckdb_sample_batch_genome_matrices

            def interrupt_third_batch(*args, **kwargs):
                if "S3" in kwargs["sample_ids"]:
                    raise RuntimeError("simulated append cancellation")
                return original_loader(*args, **kwargs)

            loader_name = "_load_duckdb_sample_batch_genome_matrices"

        monkeypatch.setattr(matrix_hdf5, loader_name, interrupt_third_batch)
        with pytest.raises(RuntimeError, match="simulated append cancellation"):
            matrix_hdf5.append_matrix_hdf5_from_duckdb(
                conn,
                sample_ids=["S2", "S3"],
                matrix_hdf5_file=matrix_file,
                export_batch_mb=0.000001,
            )

        with h5py.File(matrix_file, "r") as handle:
            assert handle.attrs[matrix_hdf5.MATRIX_BUILD_STATE_ATTR] == "appending"
            assert handle.attrs[matrix_hdf5.MATRIX_COMMITTED_SAMPLE_COUNT_ATTR] == 2
            assert handle["samples"]["sample_name"].asstr()[...].tolist() == ["S1", "S2"]
            assert handle["samples"][matrix_hdf5.MATRIX_REQUIRED_SAMPLE_DATASET].asstr()[...].tolist() == [
                "S1",
                "S2",
                "S3",
            ]
            assert handle.attrs[matrix_hdf5.MATRIX_REQUIRED_SAMPLE_COUNT_ATTR] == 3

        monkeypatch.setattr(matrix_hdf5, loader_name, original_loader)
        matrix_hdf5.append_matrix_hdf5_from_duckdb(
            conn,
            sample_ids=["S3"],
            matrix_hdf5_file=matrix_file,
            export_batch_mb=0.000001,
        )

    full = _make_database(tmp_path, name=f"append-clean-{sparse}", mode="full", cache_dir=cache_dir)
    clean_h5 = tmp_path / f"append-clean-{sparse}.h5"
    _build_matrix(
        full,
        clean_h5,
        bed_file=bed_file,
        stb_file=stb_file,
        gene_ranges=gene_ranges,
        sparse=sparse,
    )
    assert _hdf_matrix_payload(matrix_file) == _hdf_matrix_payload(clean_h5)
    with h5py.File(matrix_file, "r") as handle:
        assert handle.attrs[matrix_hdf5.MATRIX_BUILD_STATE_ATTR] == "complete"
        assert handle.attrs[matrix_hdf5.MATRIX_COMMITTED_SAMPLE_COUNT_ATTR] == len(SAMPLES)


@pytest.mark.parametrize("sparse", [False, True])
def test_matrix_append_keeps_an_eligible_sample_with_no_profile_rows_as_zero(
    tmp_path: Path,
    sparse: bool,
) -> None:
    """Stats eligibility is authoritative; an empty profile becomes an all-zero row."""
    cache_dir, bed_file, stb_file, gene_ranges = _project_files(tmp_path)
    db_file = _make_database(tmp_path, name=f"empty-row-{sparse}", mode="full", cache_dir=cache_dir)
    matrix_file = tmp_path / f"empty-row-{sparse}.h5"
    with db.connect(db_file) as conn:
        matrix_hdf5.build_matrix_hdf5_from_duckdb(
            conn,
            sample_ids=["S1", "S2"],
            output_file=matrix_file,
            genome=GENOME,
            bed_file=bed_file,
            stb_file=stb_file,
            gene_range_table=gene_ranges,
            min_cov=5,
            sparse=sparse,
        )
        conn.execute("DELETE FROM profile_positions WHERE sample_id = 'S3'")
        matrix_hdf5.append_matrix_hdf5_from_duckdb(
            conn,
            sample_ids=["S3"],
            matrix_hdf5_file=matrix_file,
        )

    with h5py.File(matrix_file, "r") as handle:
        assert handle["samples"]["sample_name"].asstr()[...].tolist() == list(SAMPLES)
        matrix = handle["matrices"]["0"]
        if isinstance(matrix, h5py.Group):
            assert int(matrix["indptr"][-1]) == int(matrix["indptr"][-2])
        else:
            assert not np.any(matrix[-1])


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


def test_group_codes_match_per_row_string_grouping() -> None:
    """Dictionary codes must find the same scaffold boundaries as per-row strings."""
    chroms = ["c1"] * 3 + ["c2"] * 2 + ["c1"] * 4
    genomes = ["g1"] * 5 + ["g2"] * 4
    batch = pa.RecordBatch.from_arrays(
        [pa.array(chroms), pa.array(genomes)], ["chrom", "genome"]
    )
    chrom_codes, chrom_values = allele_mask._group_codes(batch, "chrom")
    genome_codes, genome_values = allele_mask._group_codes(batch, "genome")

    expected = np.flatnonzero(
        (np.asarray(chroms[1:], dtype=object) != np.asarray(chroms[:-1], dtype=object))
        | (np.asarray(genomes[1:], dtype=object) != np.asarray(genomes[:-1], dtype=object))
    )
    actual = np.flatnonzero(
        (chrom_codes[1:] != chrom_codes[:-1]) | (genome_codes[1:] != genome_codes[:-1])
    )
    assert actual.tolist() == expected.tolist()
    assert [chrom_values[code] for code in chrom_codes] == chroms
    assert [genome_values[code] for code in genome_codes] == genomes


def test_group_codes_spell_nulls_as_the_skipped_genome_marker() -> None:
    """A null genome must still decode to 'None', which the import skip-list uses."""
    batch = pa.RecordBatch.from_arrays(
        [pa.array(["c1", "c1"]), pa.array([None, "g1"])], ["chrom", "genome"]
    )
    codes, values = allele_mask._group_codes(batch, "genome")
    assert values[codes[0]] == "None"
    assert values[codes[1]] == "g1"


def test_dictionary_encoded_input_groups_identically(tmp_path: Path) -> None:
    """Parquet may hand back dictionary columns; grouping must not care."""
    plain = pa.RecordBatch.from_arrays(
        [pa.array(["c1", "c1", "c2"]), pa.array(["g1"] * 3)], ["chrom", "genome"]
    )
    encoded = pa.RecordBatch.from_arrays(
        [plain.column("chrom").dictionary_encode(), plain.column("genome")],
        ["chrom", "genome"],
    )
    plain_codes, plain_values = allele_mask._group_codes(plain, "chrom")
    encoded_codes, encoded_values = allele_mask._group_codes(encoded, "chrom")
    assert [plain_values[c] for c in plain_codes] == [encoded_values[c] for c in encoded_codes]


def test_blocks_flush_in_batches_without_changing_stored_bytes(
    tmp_path: Path, monkeypatch
) -> None:
    """Forcing a flush after every block must produce identical rows."""
    cache_dir, _bed, _stb, _genes = _project_files(tmp_path)

    def stored_blocks(flush_rows: int) -> list[tuple]:
        monkeypatch.setattr(allele_mask, "BLOCK_FLUSH_MAX_ROWS", flush_rows)
        db_file = _make_database(
            tmp_path, name=f"flush{flush_rows}", mode="allele-mask", cache_dir=cache_dir
        )
        with db.connect(db_file) as conn:
            return conn.execute(
                "SELECT sample_id, segment_id, presence, deviation, covered_positions, "
                "payload_hash FROM allele_mask_profile_blocks ORDER BY sample_id, segment_id"
            ).fetchall()

    one_at_a_time = stored_blocks(1)
    batched = stored_blocks(2048)
    assert one_at_a_time == batched
    assert len(batched) > 1


def test_reference_cache_is_reused_across_samples(tmp_path: Path, monkeypatch) -> None:
    """One importer must decode each reference scaffold once, not once per sample."""
    cache_dir, _bed, _stb, _genes = _project_files(tmp_path)
    loads: list[tuple[str, str]] = []
    original = allele_mask.load_reference_segment

    def counting_load(conn, *, genome, chrom):
        loads.append((genome, chrom))
        return original(conn, genome=genome, chrom=chrom)

    monkeypatch.setattr(allele_mask, "load_reference_segment", counting_load)
    shared = allele_mask.AlleleMaskWriteCache()
    db_file = tmp_path / "shared.duckdb"
    with db.connect(db_file) as conn:
        db.configure_profile_storage(conn, mode="allele-mask", min_cov=5)
        for sample in SAMPLES:
            db.add_runs(conn, [sample])
            db.import_profile_bundle(
                conn,
                _write_bundle(tmp_path, sample),
                cache_dir=cache_dir,
                reference_cache=shared,
            )
    # Two scaffolds, three samples: without the shared cache this would be six.
    assert sorted(loads) == [(GENOME, "contigA"), (GENOME, "contigB")]


def test_reference_cache_evicts_beyond_its_base_budget(tmp_path: Path) -> None:
    """The cache must stay bounded so a wide genome set cannot grow it forever."""
    cache_dir, _bed, _stb, _genes = _project_files(tmp_path)
    db_file = _make_database(tmp_path, name="evict", mode="allele-mask", cache_dir=cache_dir)
    with db.connect(db_file) as conn:
        tiny = allele_mask.AlleleMaskWriteCache(max_cached_bases=1)
        first = tiny.segment(conn, genome=GENOME, chrom="contigA")
        tiny.segment(conn, genome=GENOME, chrom="contigB")
        assert len(tiny._segments) == 1
        again = tiny.segment(conn, genome=GENOME, chrom="contigA")
        assert np.array_equal(first.masks, again.masks)


def test_allele_mask_build_fetches_one_duckdb_block_per_matrix_batch(
    tmp_path: Path, monkeypatch
) -> None:
    """A matrix batch must be backed by one bulk DuckDB block fetch."""
    cache_dir, bed_file, stb_file, gene_ranges = _project_files(tmp_path)
    db_file = _make_database(tmp_path, name="batch-query", mode="allele-mask", cache_dir=cache_dir)
    calls: list[list[str]] = []
    original_fetch = allele_mask.fetch_profile_block_rows

    def recording_fetch(*args, **kwargs):
        calls.append(list(kwargs["sample_ids"]))
        return original_fetch(*args, **kwargs)

    monkeypatch.setattr(allele_mask, "fetch_profile_block_rows", recording_fetch)
    with db.connect(db_file) as conn:
        matrix_hdf5.build_matrix_hdf5_from_duckdb(
            conn,
            sample_ids=list(SAMPLES),
            output_file=tmp_path / "batch-query.h5",
            genome=GENOME,
            bed_file=bed_file,
            stb_file=stb_file,
            gene_range_table=gene_ranges,
            min_cov=5,
            export_batch_mb=128,
        )
    assert calls == [list(SAMPLES)]


def test_sparse_allele_mask_build_bypasses_dense_batch_materialization(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    cache_dir, bed_file, stb_file, gene_ranges = _project_files(tmp_path)
    db_file = _make_database(tmp_path, name="direct-sparse", mode="allele-mask", cache_dir=cache_dir)

    def reject_dense_loader(*_args, **_kwargs):
        raise AssertionError("sparse allele-mask build used the dense batch loader")

    monkeypatch.setattr(
        matrix_hdf5,
        "_load_duckdb_sample_batch_genome_matrices",
        reject_dense_loader,
    )
    matrix_file = tmp_path / "direct-sparse.h5"
    with db.connect(db_file) as conn:
        matrix_hdf5.build_matrix_hdf5_from_duckdb(
            conn,
            sample_ids=list(SAMPLES),
            output_file=matrix_file,
            genome=GENOME,
            bed_file=bed_file,
            stb_file=stb_file,
            gene_range_table=gene_ranges,
            min_cov=5,
            sparse=True,
        )
    assert _hdf_matrix_payload(matrix_file)[0][-1] > 0


def test_sparse_allele_mask_build_commits_bounded_stream_windows(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    cache_dir, bed_file, stb_file, gene_ranges = _project_files(tmp_path)
    db_file = _make_database(tmp_path, name="stream-checkpoints", mode="allele-mask", cache_dir=cache_dir)
    monkeypatch.setattr(matrix_hdf5, "SPARSE_STREAM_MAX_SAMPLES", 2)
    monkeypatch.setattr(matrix_hdf5, "SPARSE_STREAM_TARGET_BYTES", 1024**3)
    events: list[dict[str, object]] = []
    matrix_file = tmp_path / "stream-checkpoints.h5"
    with db.connect(db_file) as conn:
        matrix_hdf5.build_matrix_hdf5_from_duckdb(
            conn,
            sample_ids=list(SAMPLES),
            output_file=matrix_file,
            genome=GENOME,
            bed_file=bed_file,
            stb_file=stb_file,
            gene_range_table=gene_ranges,
            min_cov=5,
            sparse=True,
            progress_callback=events.append,
        )

    checkpoints = [event for event in events if event["phase"] == "checkpoint"]
    assert [event["completed"] for event in checkpoints] == [2, 3]
    assert all("query=" in str(event["detail"]) for event in checkpoints)
    assert all("decode=" in str(event["detail"]) for event in checkpoints)
    assert all("write=" in str(event["detail"]) for event in checkpoints)


def test_sparse_allele_mask_build_resumes_from_stream_checkpoint(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    cache_dir, bed_file, stb_file, gene_ranges = _project_files(tmp_path)
    db_file = _make_database(tmp_path, name="stream-resume", mode="allele-mask", cache_dir=cache_dir)
    matrix_file = tmp_path / "stream-resume.h5"
    checkpoint_file = matrix_file.with_suffix(".h5.tmp")
    monkeypatch.setattr(matrix_hdf5, "SPARSE_STREAM_MAX_SAMPLES", 1)
    original_loader = matrix_hdf5._load_allele_mask_sample_genome_sparse

    def interrupt_second_sample(*args, **kwargs):
        if kwargs["sample_id"] == "S2":
            raise RuntimeError("simulated sparse cancellation")
        return original_loader(*args, **kwargs)

    monkeypatch.setattr(
        matrix_hdf5,
        "_load_allele_mask_sample_genome_sparse",
        interrupt_second_sample,
    )
    with db.connect(db_file) as conn:
        with pytest.raises(RuntimeError, match="simulated sparse cancellation"):
            matrix_hdf5.build_matrix_hdf5_from_duckdb(
                conn,
                sample_ids=list(SAMPLES),
                output_file=matrix_file,
                genome=GENOME,
                bed_file=bed_file,
                stb_file=stb_file,
                gene_range_table=gene_ranges,
                min_cov=5,
                sparse=True,
            )
    with h5py.File(checkpoint_file, "r") as handle:
        assert handle.attrs[matrix_hdf5.MATRIX_COMMITTED_SAMPLE_COUNT_ATTR] == 1
        assert handle["samples"]["sample_name"].asstr()[...].tolist() == ["S1"]

    monkeypatch.setattr(
        matrix_hdf5,
        "_load_allele_mask_sample_genome_sparse",
        original_loader,
    )
    with db.connect(db_file) as conn:
        matrix_hdf5.build_matrix_hdf5_from_duckdb(
            conn,
            sample_ids=list(SAMPLES),
            output_file=matrix_file,
            genome=GENOME,
            bed_file=bed_file,
            stb_file=stb_file,
            gene_range_table=gene_ranges,
            min_cov=5,
            sparse=True,
        )
    clean_file = tmp_path / "stream-resume-clean.h5"
    _build_matrix(
        db_file,
        clean_file,
        bed_file=bed_file,
        stb_file=stb_file,
        gene_ranges=gene_ranges,
        sparse=True,
    )
    assert _hdf_matrix_payload(matrix_file) == _hdf_matrix_payload(clean_file)


def test_export_batch_mb_divides_samples_into_ordered_matrix_batches(
    tmp_path: Path, monkeypatch
) -> None:
    """The memory estimate must directly determine ordered DuckDB/HDF batches."""
    cache_dir, bed_file, stb_file, gene_ranges = _project_files(tmp_path)
    db_file = _make_database(tmp_path, name="ordered-batches", mode="allele-mask", cache_dir=cache_dir)
    calls: list[list[str]] = []
    original_loader = matrix_hdf5._load_duckdb_sample_batch_genome_matrices

    def recording_loader(*args, **kwargs):
        calls.append(list(kwargs["sample_ids"]))
        return original_loader(*args, **kwargs)

    monkeypatch.setattr(matrix_hdf5, "_load_duckdb_sample_batch_genome_matrices", recording_loader)
    # This reference has 16 matrix positions, so 32 bytes targets two bitmask samples.
    two_sample_batch_mb = 32 / (1024**2)
    matrix_file = tmp_path / "ordered-batches.h5"
    with db.connect(db_file) as conn:
        matrix_hdf5.build_matrix_hdf5_from_duckdb(
            conn,
            sample_ids=list(SAMPLES),
            output_file=matrix_file,
            genome=GENOME,
            bed_file=bed_file,
            stb_file=stb_file,
            gene_range_table=gene_ranges,
            min_cov=5,
            export_batch_mb=two_sample_batch_mb,
        )
    assert calls == [["S1", "S2"], ["S3"]]
    with h5py.File(matrix_file, "r") as handle:
        assert handle["samples"]["sample_name"].asstr()[...].tolist() == list(SAMPLES)


def test_fetch_profile_block_rows_groups_every_requested_sample(tmp_path: Path) -> None:
    """Samples with no blocks must still appear, so callers can tell empty from absent."""
    cache_dir, _bed, _stb, _genes = _project_files(tmp_path)
    db_file = _make_database(tmp_path, name="fetch", mode="allele-mask", cache_dir=cache_dir)
    with db.connect(db_file) as conn:
        segments = [
            int(row[0])
            for row in conn.execute(
                "SELECT segment_id FROM allele_mask_reference_segments"
            ).fetchall()
        ]
        grouped = allele_mask.fetch_profile_block_rows(
            conn, sample_ids=["S1", "S2", "missing"], segment_ids=segments
        )
    assert set(grouped) == {"S1", "S2", "missing"}
    assert grouped["missing"] == []
    assert grouped["S1"]


def test_fetch_profile_block_rows_adds_scalar_bounds_for_duckdb_pruning() -> None:
    calls: list[tuple[str, list[object]]] = []

    class Result:
        def fetchmany(self, _size):
            return []

    class Connection:
        def execute(self, query, parameters):
            calls.append((query, parameters))
            return Result()

    grouped = allele_mask.fetch_profile_block_rows(
        Connection(),
        sample_ids=["S9", "S1"],
        segment_ids=[9, 2],
    )

    query, parameters = calls[0]
    assert "sample_id BETWEEN ? AND ?" in query
    assert "segment_id BETWEEN ? AND ?" in query
    assert parameters[:4] == ["S1", "S9", 2, 9]
    assert set(grouped) == {"S1", "S9"}


def test_profile_block_stream_uses_bounded_fetches_and_yields_empty_samples() -> None:
    rows = [
        ("S1", 2, b"p1", b"d1", 1, "h1"),
        ("S3", 2, b"p3", b"d3", 1, "h3"),
    ]

    class Result:
        def __init__(self):
            self.offset = 0
            self.fetch_sizes: list[int] = []

        def fetchmany(self, size):
            self.fetch_sizes.append(size)
            chunk = rows[self.offset : self.offset + 1]
            self.offset += len(chunk)
            return chunk

    result = Result()

    class Connection:
        def execute(self, _query, _parameters):
            return result

    metrics: dict[str, float | int] = {}
    streamed = list(
        allele_mask.iter_profile_block_rows(
            Connection(),
            sample_ids=["S3", "S1", "S2"],
            segment_ids=[2],
            fetch_rows=1,
            metrics=metrics,
        )
    )
    assert [sample_id for sample_id, _rows in streamed] == ["S1", "S2", "S3"]
    assert streamed[1][1] == []
    assert result.fetch_sizes == [1, 1, 1]
    assert metrics["fetched_rows"] == 2


def test_matrix_processing_log_reports_batch_operation(capsys) -> None:
    logger = ThrottledMatrixLogger("MATRIX-BUILD")
    logger(
        {
            "phase": "processing",
            "completed": 0,
            "total": 10,
            "sample_name": "S1",
            "genome": GENOME,
            "stored_rows": 0,
            "detail": "fetching batch_samples=4 first_sample=S1",
        }
    )
    captured = capsys.readouterr()
    assert "MATRIX-BUILD PROCESSING" in captured.err
    assert "detail=fetching batch_samples=4 first_sample=S1" in captured.err


def test_full_profile_matrix_build_reports_loader_stage_timings(tmp_path: Path) -> None:
    cache_dir, bed_file, stb_file, gene_ranges = _project_files(tmp_path)
    db_file = _make_database(tmp_path, name="full-stage-timing", mode="full", cache_dir=cache_dir)
    events: list[dict[str, object]] = []
    with db.connect(db_file) as conn:
        matrix_hdf5.build_matrix_hdf5_from_duckdb(
            conn,
            sample_ids=list(SAMPLES),
            output_file=tmp_path / "full-stage-timing.h5",
            genome=GENOME,
            bed_file=bed_file,
            stb_file=stb_file,
            gene_range_table=gene_ranges,
            min_cov=5,
            sparse=True,
            progress_callback=events.append,
        )

    details = [
        str(event["detail"])
        for event in events
        if event["phase"] == "processing" and str(event.get("detail", "")).startswith("full-profile")
    ]
    assert [detail.split("stage=", 1)[1].split(" ", 1)[0] for detail in details] == [
        "allocate",
        "query",
        "stream",
    ]
    assert all("seconds=" in detail and "rows=" in detail for detail in details)


@pytest.mark.parametrize("storage_mode", ["bitmask", "counts"])
@pytest.mark.parametrize("sparse", [False, True])
def test_full_profile_arrow_stream_matches_frame_conversion(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    storage_mode: str,
    sparse: bool,
) -> None:
    cache_dir, bed_file, stb_file, gene_ranges = _project_files(tmp_path)
    db_file = _make_database(tmp_path, name=f"arrow-{storage_mode}-{sparse}", mode="full", cache_dir=cache_dir)

    def reject_frame_conversion(*_args, **_kwargs):
        raise AssertionError("full-profile matrix build materialized a Polars sample frame")

    monkeypatch.setattr(matrix_hdf5, "_profile_frame_to_matrix", reject_frame_conversion)
    matrix_file = tmp_path / f"arrow-{storage_mode}-{sparse}.h5"
    with db.connect(db_file) as conn:
        matrix_hdf5.build_matrix_hdf5_from_duckdb(
            conn,
            sample_ids=list(SAMPLES),
            output_file=matrix_file,
            genome=GENOME,
            bed_file=bed_file,
            stb_file=stb_file,
            gene_range_table=gene_ranges,
            storage_mode=storage_mode,
            count_dtype="uint16" if storage_mode == "counts" else None,
            min_cov=5,
            sparse=sparse,
        )

    with h5py.File(matrix_file, "r") as handle:
        sample_names = handle["samples"]["sample_name"].asstr()[...].tolist()
        node = handle["matrices"]["0"]
        if isinstance(node, h5py.Group):
            dense = np.zeros(
                (len(SAMPLES), 16) if storage_mode == "bitmask" else (len(SAMPLES), 16, 4),
                dtype=np.uint8 if storage_mode == "bitmask" else np.uint16,
            )
            indptr = node["indptr"][...]
            indices = node["indices"][...]
            values = node["values"][...]
            flat = dense.reshape((len(SAMPLES), -1))
            for row_index in range(len(SAMPLES)):
                start, stop = int(indptr[row_index]), int(indptr[row_index + 1])
                flat[row_index, indices[start:stop]] = values[start:stop]
        else:
            dense = node[...]

    expected = np.zeros_like(dense)
    for row_index, sample_id in enumerate(SAMPLES):
        for chrom, pos, a_count, c_count, g_count, t_count in _profile_rows(sample_id):
            if a_count + c_count + g_count + t_count < 5:
                continue
            axis = pos - 1 if chrom == "contigA" else 11 + pos - 1
            counts = np.asarray([a_count, t_count, c_count, g_count], dtype=np.uint16)
            if storage_mode == "bitmask":
                expected[row_index, axis] = int(((counts > 0) * matrix_hdf5.BITMASK_BASE_BITS).sum())
            else:
                expected[row_index, axis, :] = counts
    assert sample_names == list(SAMPLES)
    assert np.array_equal(dense, expected)


def test_full_profile_arrow_stream_crosses_record_batch_boundary(tmp_path: Path) -> None:
    genome = "GCF_ARROW.1"
    length = 70_000
    db_file = tmp_path / "arrow-boundary.duckdb"
    with db.connect(db_file) as conn:
        conn.execute(
            f"""
            INSERT INTO profile_positions
            SELECT 'S1', 'contig', pos, '{genome}',
                   CAST(5 AS USMALLINT), CAST(0 AS USMALLINT),
                   CAST(0 AS USMALLINT), CAST(0 AS USMALLINT), NULL
            FROM range(1, {length + 1}) positions(pos)
            """
        )
    bed_file = tmp_path / "arrow-boundary.bed"
    bed_file.write_text(f"contig\t0\t{length}\n")
    stb_file = tmp_path / "arrow-boundary.stb"
    stb_file.write_text(f"contig\t{genome}\n")
    matrix_file = tmp_path / "arrow-boundary.h5"

    with db.connect(db_file) as conn:
        matrix_hdf5.build_matrix_hdf5_from_duckdb(
            conn,
            sample_ids=["S1"],
            output_file=matrix_file,
            genome=genome,
            bed_file=bed_file,
            stb_file=stb_file,
            min_cov=5,
            sparse=False,
        )
    with h5py.File(matrix_file, "r") as handle:
        matrix = handle["matrices"]["0"][0]
        assert matrix.shape == (length,)
        assert np.all(matrix == 1)


def test_publish_moves_outputs_instead_of_copying(tmp_path: Path) -> None:
    """Publishing renames, so the profile parquet is not written twice."""
    source = tmp_path / "work" / "out.parquet"
    source.parent.mkdir()
    source.write_bytes(b"payload")
    destination = tmp_path / "published" / "out.parquet"

    workflows._atomic_publish(source, destination)

    assert destination.read_bytes() == b"payload"
    assert not source.exists()


def test_publish_falls_back_to_copying_across_filesystems(tmp_path: Path, monkeypatch) -> None:
    """Node-local scratch and networked output are different devices."""
    source = tmp_path / "work" / "out.parquet"
    source.parent.mkdir()
    source.write_bytes(b"payload")
    destination = tmp_path / "published" / "out.parquet"

    real_replace = workflows.os.replace
    calls = {"count": 0}

    def cross_device(*args, **kwargs):
        calls["count"] += 1
        if calls["count"] == 1:
            raise OSError(errno.EXDEV, "Cross-device link")
        return real_replace(*args, **kwargs)

    monkeypatch.setattr(workflows.os, "replace", cross_device)
    workflows._atomic_publish(source, destination)

    assert destination.read_bytes() == b"payload"
    assert source.exists()  # the copy path leaves the caller to clean up


def test_publish_refuses_an_empty_output(tmp_path: Path) -> None:
    """An empty output must never be published as if it were complete."""
    source = tmp_path / "empty.parquet"
    source.write_bytes(b"")
    with pytest.raises(RuntimeError, match="empty or missing"):
        workflows._atomic_publish(source, tmp_path / "published.parquet")


def test_plan_uses_only_stats_and_hdf5_membership(tmp_path: Path) -> None:
    """Registry rows must not influence whether a genome needs matrix work."""
    cache_dir, bed_file, stb_file, gene_ranges = _project_files(tmp_path)
    db_file = _make_database(tmp_path, name="plan-hdf", mode="allele-mask", cache_dir=cache_dir)
    matrix_dir = tmp_path / "matrices"
    matrix_dir.mkdir()
    _build_matrix(
        db_file,
        matrix_dir / f"{GENOME}.h5",
        bed_file=bed_file,
        stb_file=stb_file,
        gene_ranges=gene_ranges,
        sparse=False,
    )
    with db.connect(db_file) as conn:
        conn.execute("DELETE FROM matrix_store_samples")
        conn.execute("DELETE FROM matrix_stores")
        plan = workflows.plan_matrix_sync_build(
            conn,
            matrix_dir=matrix_dir,
            genomes=[GENOME],
            filters=db.MatrixFilters(),
        )
    assert plan.pending == []
    assert plan.up_to_date == (GENOME,)


def test_plan_marks_a_stats_eligible_sample_missing_from_hdf_as_behind(tmp_path: Path) -> None:
    """Eligibility comes from stats; membership comes only from the HDF sample axis."""
    cache_dir, bed_file, stb_file, gene_ranges = _project_files(tmp_path)
    db_file = _make_database(tmp_path, name="plan-behind", mode="allele-mask", cache_dir=cache_dir)
    matrix_dir = tmp_path / "matrices"
    matrix_dir.mkdir()
    matrix_file = matrix_dir / f"{GENOME}.h5"
    with db.connect(db_file) as conn:
        matrix_hdf5.build_matrix_hdf5_from_duckdb(
            conn,
            sample_ids=["S1", "S2"],
            output_file=matrix_file,
            genome=GENOME,
            bed_file=bed_file,
            stb_file=stb_file,
            gene_range_table=gene_ranges,
            min_cov=5,
        )
        plan = workflows.plan_matrix_sync_build(
            conn,
            matrix_dir=matrix_dir,
            genomes=[GENOME],
            filters=db.MatrixFilters(),
        )
    assert plan.behind == (GENOME,)
    assert plan.pending == [GENOME]


def test_plan_does_not_scan_profiles_to_validate_eligible_samples(tmp_path: Path) -> None:
    """Missing profile rows are represented by zeros rather than a terabyte-table scan."""
    cache_dir, bed_file, stb_file, gene_ranges = _project_files(tmp_path)
    db_file = _make_database(tmp_path, name="plan-no-profile-scan", mode="full", cache_dir=cache_dir)
    matrix_dir = tmp_path / "matrices"
    matrix_dir.mkdir()
    matrix_file = matrix_dir / f"{GENOME}.h5"
    with db.connect(db_file) as conn:
        matrix_hdf5.build_matrix_hdf5_from_duckdb(
            conn,
            sample_ids=["S1", "S2"],
            output_file=matrix_file,
            genome=GENOME,
            bed_file=bed_file,
            stb_file=stb_file,
            gene_range_table=gene_ranges,
            min_cov=5,
        )
        conn.execute("DELETE FROM profile_positions WHERE sample_id = 'S3'")
        plan = workflows.plan_matrix_sync_build(
            conn,
            matrix_dir=matrix_dir,
            genomes=[GENOME],
            filters=db.MatrixFilters(),
        )
    assert plan.behind == (GENOME,)


def test_plan_skips_genomes_without_stats_eligible_samples(tmp_path: Path) -> None:
    cache_dir, _bed_file, _stb_file, _gene_ranges = _project_files(tmp_path)
    db_file = _make_database(tmp_path, name="plan-filtered", mode="allele-mask", cache_dir=cache_dir)
    with db.connect(db_file) as conn:
        plan = workflows.plan_matrix_sync_build(
            conn,
            matrix_dir=tmp_path / "matrices",
            genomes=[GENOME],
            filters=db.MatrixFilters(min_coverage=100.0),
        )
    assert plan.skipped == (GENOME,)
    assert plan.pending == []


def test_plan_marks_an_eligible_genome_without_a_file_as_missing(tmp_path: Path) -> None:
    cache_dir, _bed_file, _stb_file, _gene_ranges = _project_files(tmp_path)
    db_file = _make_database(tmp_path, name="plan-missing", mode="allele-mask", cache_dir=cache_dir)
    with db.connect(db_file) as conn:
        plan = workflows.plan_matrix_sync_build(
            conn,
            matrix_dir=tmp_path / "matrices",
            genomes=[GENOME],
            filters=db.MatrixFilters(),
        )
    assert plan.missing == (GENOME,)
    assert plan.pending == [GENOME]


def _strip_checkpoint_attrs(matrix_file: Path) -> None:
    with h5py.File(matrix_file, "r+") as handle:
        del handle.attrs[matrix_hdf5.MATRIX_BUILD_STATE_ATTR]
        del handle.attrs[matrix_hdf5.MATRIX_COMMITTED_SAMPLE_COUNT_ATTR]
        del handle.attrs[matrix_hdf5.MATRIX_REQUIRED_SAMPLE_COUNT_ATTR]
        del handle.attrs[matrix_hdf5.MATRIX_REQUIRED_SAMPLE_DIGEST_ATTR]
        del handle["samples"][matrix_hdf5.MATRIX_REQUIRED_SAMPLE_DATASET]


def test_plan_accepts_and_upgrades_legacy_hdf_checkpoint_metadata(tmp_path: Path) -> None:
    """Old matrices remain usable and pay the full sample-list read only once."""
    cache_dir, bed_file, stb_file, gene_ranges = _project_files(tmp_path)
    db_file = _make_database(tmp_path, name="plan-legacy", mode="allele-mask", cache_dir=cache_dir)
    matrix_dir = tmp_path / "matrices"
    matrix_dir.mkdir()
    matrix_file = matrix_dir / f"{GENOME}.h5"
    _build_matrix(
        db_file,
        matrix_file,
        bed_file=bed_file,
        stb_file=stb_file,
        gene_ranges=gene_ranges,
        sparse=False,
    )
    _strip_checkpoint_attrs(matrix_file)

    with db.connect(db_file) as conn:
        plan = workflows.plan_matrix_sync_build(
            conn,
            matrix_dir=matrix_dir,
            genomes=[GENOME],
            filters=db.MatrixFilters(),
        )
    assert plan.up_to_date == (GENOME,)
    with h5py.File(matrix_file, "r") as handle:
        assert handle.attrs[matrix_hdf5.MATRIX_BUILD_STATE_ATTR] == "complete"
        assert handle.attrs[matrix_hdf5.MATRIX_COMMITTED_SAMPLE_COUNT_ATTR] == len(SAMPLES)
        assert handle.attrs[matrix_hdf5.MATRIX_REQUIRED_SAMPLE_COUNT_ATTR] == len(SAMPLES)
        assert handle.attrs[matrix_hdf5.MATRIX_REQUIRED_SAMPLE_DIGEST_ATTR] == matrix_hdf5.sample_id_set_digest(
            list(SAMPLES)
        )
        assert handle["samples"][matrix_hdf5.MATRIX_REQUIRED_SAMPLE_DATASET].asstr()[...].tolist() == sorted(SAMPLES)


def test_eligibility_summary_digest_matches_hdf_required_sample_digest(tmp_path: Path) -> None:
    cache_dir, bed_file, stb_file, gene_ranges = _project_files(tmp_path)
    db_file = _make_database(tmp_path, name="summary-digest", mode="allele-mask", cache_dir=cache_dir)
    matrix_file = tmp_path / f"{GENOME}.h5"
    _build_matrix(
        db_file,
        matrix_file,
        bed_file=bed_file,
        stb_file=stb_file,
        gene_ranges=gene_ranges,
        sparse=True,
    )

    with db.connect(db_file) as conn:
        summaries = db.eligible_genome_summaries(
            conn,
            genomes=[GENOME],
            filters=db.MatrixFilters(),
        )
    assert len(summaries) == 1
    with h5py.File(matrix_file, "r") as handle:
        assert summaries[0].sample_count == handle.attrs[matrix_hdf5.MATRIX_REQUIRED_SAMPLE_COUNT_ATTR]
        assert summaries[0].sample_digest == handle.attrs[matrix_hdf5.MATRIX_REQUIRED_SAMPLE_DIGEST_ATTR]


def test_append_rejects_filters_that_remove_committed_matrix_samples(tmp_path: Path) -> None:
    cache_dir, bed_file, stb_file, gene_ranges = _project_files(tmp_path)
    db_file = _make_database(tmp_path, name="removed-required", mode="allele-mask", cache_dir=cache_dir)
    matrix_file = tmp_path / f"{GENOME}.h5"
    _build_matrix(
        db_file,
        matrix_file,
        bed_file=bed_file,
        stb_file=stb_file,
        gene_ranges=gene_ranges,
        sparse=True,
    )

    with db.connect(db_file) as conn:
        with pytest.raises(ValueError, match="filters remove samples already stored"):
            workflows.append_matrix_from_database(
                conn,
                matrix_file=matrix_file,
                filters=db.MatrixFilters(min_coverage=1000),
                register=False,
            )

from __future__ import annotations

import gzip
import json
from pathlib import Path

from click.testing import CliRunner
import duckdb
import h5py
import numpy as np
import polars as pl

from metatrawl import cli
from metatrawl import db
from metatrawl import genome_views


def _project_database(path: Path) -> Path:
    with db.connect(path) as conn:
        now = 1.0
        for sample, coverage, breadth, ber, ref_ani, abundance in (
            ("sample_a", 3.0, 0.91, 0.81, 0.999, 0.15),
            ("sample_b", 2.0, 0.82, 0.77, 0.998, 0.12),
            ("sample_c", 4.0, 0.95, 0.88, 0.997, 0.20),
        ):
            conn.execute("INSERT INTO sra_runs VALUES (?, 'complete', ?, ?, NULL)", [sample, now, now])
            conn.execute("INSERT INTO samples VALUES (?, ?, 'complete', ?, ?)", [sample, sample, now, now])
            conn.execute(
                """INSERT INTO genome_stats
                   (sample_id, genome, coverage, breadth, ber, ref_ani)
                   VALUES (?, 'genome_a', ?, ?, ?, ?)""",
                [sample, coverage, breadth, ber, ref_ani],
            )
            conn.execute(
                "INSERT INTO sylph_abundance VALUES (?, 'genome_a', 'genome_a', ?)",
                [sample, abundance],
            )
    return path


def _compare_database(
    path: Path,
    *,
    completed_pairs: int = 3,
    legacy_ani: bool = False,
    include_catalogs: bool = True,
    include_checkpoint: bool = True,
    result_pairs: int = 3,
) -> Path:
    with duckdb.connect(str(path)) as conn:
        if include_catalogs:
            conn.execute("CREATE TABLE matrix_compare_samples (sample_idx INTEGER, sample_name VARCHAR)")
            conn.execute("CREATE TABLE matrix_compare_genomes (genome_idx INTEGER, genome VARCHAR)")
        if include_checkpoint:
            conn.execute(
                "CREATE TABLE matrix_compare_completed_pair_genomes "
                "(sample_idx_1 INTEGER, sample_idx_2 INTEGER, genome_idx INTEGER)"
            )
        ani_column = "genome_pop_ani" if legacy_ani else "genome_ani"
        conn.execute(
            f"""
            CREATE TABLE matrix_compare_results (
              sample_idx_1 INTEGER,
              sample_idx_2 INTEGER,
              sample_1 VARCHAR,
              sample_2 VARCHAR,
              genome_idx INTEGER,
              genome VARCHAR,
              total_positions BIGINT,
              share_allele_pos BIGINT,
              {ani_column} DOUBLE,
              max_consecutive_length BIGINT
            )
            """
        )
        if include_catalogs:
            conn.executemany(
                "INSERT INTO matrix_compare_samples VALUES (?, ?)",
                [(0, "sample_a"), (1, "sample_b"), (2, "sample_c")],
            )
            conn.execute("INSERT INTO matrix_compare_genomes VALUES (0, 'genome_a')")
        if include_checkpoint:
            completed = [(0, 1, 0), (0, 2, 0), (1, 2, 0)][:completed_pairs]
            conn.executemany("INSERT INTO matrix_compare_completed_pair_genomes VALUES (?, ?, ?)", completed)
        result_rows = [
            (0, 1, "sample_a", "sample_b", 120_000, 119_940, 99.95, 11_000),
            (0, 2, "sample_a", "sample_c", 110_000, 109_780, 99.80, 8_000),
            (1, 2, "sample_b", "sample_c", 100_000, 99_900, 99.90, 9_000),
        ]
        conn.executemany(
            "INSERT INTO matrix_compare_results VALUES (?, ?, ?, ?, 0, 'genome_a', ?, ?, ?, ?)",
            result_rows[:result_pairs],
        )
    return path


def _fraction_compare_database(path: Path) -> Path:
    samples = [f"sample_{letter}" for letter in "abcdef"]
    qualifying_edges = [
        (0, 1),
        (0, 2),
        (0, 3),
        (1, 2),
        (1, 3),
        (2, 3),
        (0, 4),
        (1, 5),
    ]
    with duckdb.connect(str(path)) as conn:
        conn.execute("CREATE TABLE matrix_compare_samples (sample_idx INTEGER, sample_name VARCHAR)")
        conn.execute("CREATE TABLE matrix_compare_genomes (genome_idx INTEGER, genome VARCHAR)")
        conn.execute(
            "CREATE TABLE matrix_compare_completed_pair_genomes "
            "(sample_idx_1 INTEGER, sample_idx_2 INTEGER, genome_idx INTEGER)"
        )
        conn.execute(
            """
            CREATE TABLE matrix_compare_results (
              sample_idx_1 INTEGER,
              sample_idx_2 INTEGER,
              sample_1 VARCHAR,
              sample_2 VARCHAR,
              genome_idx INTEGER,
              genome VARCHAR,
              total_positions BIGINT,
              share_allele_pos BIGINT,
              genome_ani DOUBLE,
              max_consecutive_length BIGINT
            )
            """
        )
        conn.executemany(
            "INSERT INTO matrix_compare_samples VALUES (?, ?)",
            list(enumerate(samples)),
        )
        conn.execute("INSERT INTO matrix_compare_genomes VALUES (0, 'genome_a')")
        all_edges = [(left, right, 0) for left in range(6) for right in range(left + 1, 6)]
        conn.executemany(
            "INSERT INTO matrix_compare_completed_pair_genomes VALUES (?, ?, ?)",
            all_edges,
        )
        conn.executemany(
            "INSERT INTO matrix_compare_results VALUES (?, ?, ?, ?, 0, 'genome_a', 20000, 19980, 99.9, 10000)",
            [(left, right, samples[left], samples[right]) for left, right in qualifying_edges],
        )
    return path


def test_sync_genome_views_writes_self_contained_bundle_and_resumes(tmp_path: Path) -> None:
    project_db = _project_database(tmp_path / "metatrawl.duckdb")
    compare_dir = tmp_path / "compares"
    compare_dir.mkdir()
    _compare_database(compare_dir / "genome_a.duckdb")
    view_dir = tmp_path / "views"

    summary = genome_views.sync_genome_views(
        db_file=project_db,
        compare_dir=compare_dir,
        view_dir=view_dir,
        options=genome_views.GenomeViewOptions(neighbor_k=1),
    )

    assert summary.generated == 1
    assert summary.failed == 0
    catalog = json.loads((view_dir / "catalog.json").read_text())
    assert catalog["genome_count"] == 1
    assert catalog["genomes"][0]["path"] == "genome_a"
    assert catalog["genomes"][0]["manifest"] == "genome_a/manifest.json"
    bundle = view_dir / "genome_a"
    assert {path.name for path in bundle.iterdir()} == set(genome_views.VIEW_ARTIFACTS)
    manifest = json.loads((bundle / "manifest.json").read_text())
    assert manifest["schema_version"] == 2
    assert manifest["genome"] == "genome_a"
    assert manifest["sample_count"] == 3
    assert manifest["source_signature"]["completed_pairs"] == 3
    assert manifest["files"]["similarity_matrix"]["shape"] == [3]
    assert manifest["files"]["similarity_matrix"]["matrix_shape"] == [3, 3]
    assert manifest["files"]["similarity_matrix"]["dtype"] == "float32"
    assert manifest["files"]["total_positions_matrix"]["dtype"] == "uint64"
    for descriptor in manifest["files"].values():
        assert descriptor["size_bytes"] == (bundle / descriptor["path"]).stat().st_size

    with h5py.File(bundle / "view_data.h5", "r") as handle:
        assert handle["similarity_ani"].shape == (3, 3)
        assert handle["total_positions"].shape == (3, 3)
        assert handle["linkage"].shape == (2, 4)
        hdf5_similarity = handle["similarity_ani"][...]
        hdf5_total_positions = handle["total_positions"][...]
        assert sorted(value.decode() if isinstance(value, bytes) else str(value) for value in handle["samples"][:]) == [
            "sample_a",
            "sample_b",
            "sample_c",
        ]
    from scipy.spatial.distance import squareform

    condensed_similarity = np.frombuffer(
        gzip.decompress((bundle / "similarity_ani.condensed.f32.gz").read_bytes()),
        dtype="<f4",
    )
    condensed_total_positions = np.frombuffer(
        gzip.decompress((bundle / "total_positions.condensed.u64.gz").read_bytes()),
        dtype="<u8",
    )
    web_similarity = squareform(condensed_similarity)
    np.fill_diagonal(web_similarity, 100.0)
    web_total_positions = squareform(condensed_total_positions)
    np.testing.assert_array_equal(web_similarity, hdf5_similarity)
    np.testing.assert_array_equal(web_total_positions, hdf5_total_positions)
    assert web_total_positions[0, 1] == 120_000
    assert web_total_positions[1, 0] == 120_000

    samples = json.loads((bundle / "samples.json").read_text())
    assert samples["sample_count"] == 3
    assert sorted(row["sample_index"] for row in samples["samples"]) == [0, 1, 2]
    assert sorted(samples["leaf_order"]) == [0, 1, 2]

    stats = pl.read_parquet(bundle / "sample_stats.parquet")
    assert stats.height == 3
    assert stats.get_column("coverage").drop_nulls().len() == 3
    web_stats = json.loads((bundle / "sample_stats.json").read_text())
    assert web_stats["orientation"] == "columnar"
    assert web_stats["row_count"] == 3
    assert web_stats["columns"]["sample_id"] == stats.get_column("sample_id").to_list()
    assert web_stats["columns"]["coverage"] == stats.get_column("coverage").to_list()

    clusters = json.loads((bundle / "clusters.json").read_text())
    assert len(clusters["assignments"]) == 3
    assert sum(group["sample_count"] for group in clusters["clusters"]["strain"]) == 3
    assert sum(group["sample_count"] for group in clusters["clusters"]["clonal"]) == 3

    dendrogram = json.loads((bundle / "dendrogram.json").read_text())
    assert dendrogram["linkage"]["format"] == "scipy"
    assert len(dendrogram["tree"]["merges"]) == 2
    assert dendrogram["tree"]["root_id"] == 4

    network = json.loads((bundle / "neighbor_network.json").read_text())
    assert len(network["nodes"]) == 3
    assert 1 <= len(network["edges"]) <= 3
    assert all(edge["total_positions"] > 0 for edge in network["edges"])

    resumed = genome_views.sync_genome_views(
        db_file=project_db,
        compare_dir=compare_dir,
        view_dir=view_dir,
        options=genome_views.GenomeViewOptions(neighbor_k=1),
    )
    assert resumed.generated == 0
    assert resumed.up_to_date == 1


def test_sync_genome_views_refreshes_when_project_stats_change(tmp_path: Path) -> None:
    project_db = _project_database(tmp_path / "metatrawl.duckdb")
    compare_dir = tmp_path / "compares"
    compare_dir.mkdir()
    _compare_database(compare_dir / "genome_a.duckdb")
    view_dir = tmp_path / "views"

    first = genome_views.sync_genome_views(
        db_file=project_db,
        compare_dir=compare_dir,
        view_dir=view_dir,
    )
    assert first.generated == 1

    with db.connect(project_db) as conn:
        conn.execute(
            "UPDATE genome_stats SET coverage = 9.5 "
            "WHERE sample_id = 'sample_a' AND genome = 'genome_a'"
        )

    refreshed = genome_views.sync_genome_views(
        db_file=project_db,
        compare_dir=compare_dir,
        view_dir=view_dir,
    )
    assert refreshed.generated == 1
    assert refreshed.up_to_date == 0
    stats = json.loads((view_dir / "genome_a" / "sample_stats.json").read_text())
    sample_index = stats["columns"]["sample_id"].index("sample_a")
    assert stats["columns"]["coverage"][sample_index] == 9.5


def test_sync_genome_views_reads_pre_v1_genome_pop_ani_table(tmp_path: Path) -> None:
    project_db = _project_database(tmp_path / "metatrawl.duckdb")
    compare_dir = tmp_path / "compares"
    compare_dir.mkdir()
    _compare_database(compare_dir / "genome_a.duckdb", legacy_ani=True)
    view_dir = tmp_path / "views"

    summary = genome_views.sync_genome_views(
        db_file=project_db,
        compare_dir=compare_dir,
        view_dir=view_dir,
    )

    assert summary.generated == 1
    assert summary.failed == 0
    manifest = json.loads((view_dir / "genome_a" / "manifest.json").read_text())
    compare_schema = manifest["source_signature"]["comparison_schema"]
    assert compare_schema["ani_column"] == "genome_pop_ani"
    assert compare_schema["completion_source"] == "matrix_compare_completed_pair_genomes"
    matrix = np.frombuffer(
        gzip.decompress(
            (view_dir / "genome_a" / "similarity_ani.condensed.f32.gz").read_bytes()
        ),
        dtype="<f4",
    )
    np.testing.assert_allclose(matrix, [99.95, 99.8, 99.9], rtol=0, atol=1e-4)


def test_sync_genome_views_infers_legacy_catalogs_and_completion_from_results(tmp_path: Path) -> None:
    project_db = _project_database(tmp_path / "metatrawl.duckdb")
    compare_dir = tmp_path / "compares"
    compare_dir.mkdir()
    _compare_database(
        compare_dir / "genome_a.duckdb",
        legacy_ani=True,
        include_catalogs=False,
        include_checkpoint=False,
    )
    view_dir = tmp_path / "views"

    summary = genome_views.sync_genome_views(
        db_file=project_db,
        compare_dir=compare_dir,
        view_dir=view_dir,
    )

    assert summary.generated == 1
    manifest = json.loads((view_dir / "genome_a" / "manifest.json").read_text())
    compare_schema = manifest["source_signature"]["comparison_schema"]
    assert compare_schema == {
        "ani_column": "genome_pop_ani",
        "sample_catalog": False,
        "genome_catalog": False,
        "completed_pairs_table": False,
        "completion_source": "distinct_result_rows",
    }


def test_sync_genome_views_rejects_incomplete_legacy_results_without_checkpoint(tmp_path: Path) -> None:
    project_db = _project_database(tmp_path / "metatrawl.duckdb")
    compare_dir = tmp_path / "compares"
    compare_dir.mkdir()
    _compare_database(
        compare_dir / "genome_a.duckdb",
        legacy_ani=True,
        include_checkpoint=False,
        result_pairs=2,
    )

    summary = genome_views.sync_genome_views(
        db_file=project_db,
        compare_dir=compare_dir,
        view_dir=tmp_path / "views",
    )

    assert summary.skipped == 1
    assert summary.ready == 0


def test_sync_genome_views_skips_incomplete_compare_database(tmp_path: Path) -> None:
    project_db = _project_database(tmp_path / "metatrawl.duckdb")
    compare_dir = tmp_path / "compares"
    compare_dir.mkdir()
    _compare_database(compare_dir / "genome_a.duckdb", completed_pairs=2)

    summary = genome_views.sync_genome_views(
        db_file=project_db,
        compare_dir=compare_dir,
        view_dir=tmp_path / "views",
    )

    assert summary.skipped == 1
    assert summary.ready == 0
    assert summary.failed == 0
    assert not (tmp_path / "views" / "genome_a" / "manifest.json").exists()


def test_genome_view_clustering_matches_zipstrain_visualization_semantics(tmp_path: Path) -> None:
    from zipstrain import visualize

    compare_db = _compare_database(tmp_path / "genome_a.duckdb")
    with duckdb.connect(str(compare_db), read_only=True) as conn:
        prepared = genome_views._prepare_genome_view(
            conn,
            genome="genome_a",
            options=genome_views.GenomeViewOptions(),
        )
        comparisons = conn.execute(
            "SELECT sample_1, sample_2, genome, total_positions, genome_ani FROM matrix_compare_results"
        ).pl()

    expected = visualize._prepare_similarity_matrix(
        comparisons.lazy(),
        genome="genome_a",
        min_comp_len=10_000,
        impute_method=97.0,
        max_null_samples=500,
        linkage_method="average",
    )

    assert prepared.samples == expected.samples
    np.testing.assert_allclose(prepared.similarity_matrix, expected.similarity_matrix, rtol=0, atol=1e-5)
    np.testing.assert_allclose(prepared.linkage_matrix, expected.linkage_matrix, rtol=0, atol=1e-5)


def test_genome_view_filters_by_genome_specific_null_fraction(tmp_path: Path) -> None:
    compare_db = _fraction_compare_database(tmp_path / "genome_a.duckdb")
    with duckdb.connect(str(compare_db), read_only=True) as conn:
        prepared = genome_views._prepare_genome_view(
            conn,
            genome="genome_a",
            options=genome_views.GenomeViewOptions(max_null_fraction=0.4),
        )
        compact_columns = [
            row[0]
            for row in conn.execute("DESCRIBE metatrawl_view_pairs").fetchall()
        ]

    assert prepared.samples == ["sample_a", "sample_b", "sample_c", "sample_d"]
    assert prepared.similarity_matrix.shape == (4, 4)
    assert compact_columns == ["sample_idx_1", "sample_idx_2", "genome_ani", "total_positions"]


def test_genome_view_absolute_null_limit_overrides_fraction(tmp_path: Path) -> None:
    compare_db = _fraction_compare_database(tmp_path / "genome_a.duckdb")
    with duckdb.connect(str(compare_db), read_only=True) as conn:
        prepared = genome_views._prepare_genome_view(
            conn,
            genome="genome_a",
            options=genome_views.GenomeViewOptions(
                max_null_fraction=1.0,
                max_null_samples=1,
            ),
        )

    assert prepared.samples == ["sample_a", "sample_b"]


def test_sync_genome_views_cli_dispatches_one_job_per_genome(tmp_path: Path, monkeypatch) -> None:
    compare_dir = tmp_path / "compares"
    compare_dir.mkdir()
    (compare_dir / "genome_a.duckdb").touch()
    (compare_dir / "genome_b.duckdb").touch()
    workflow_config = tmp_path / "workflow.toml"
    workflow_config.write_text(
        '[stages.genome_view]\n'
        'workers = 2\n'
        'threads = 4\n'
        'execution = "slurm"\n'
        '[genome_view]\n'
        'max_null_fraction = 0.35\n'
    )
    commands: list[tuple[str, list[str], str]] = []

    class FakeRuntime:
        def __init__(self, config, *, state_dir, logger, runner=None):
            pass

        def run(self, stage, command, *, sample, stdout_file=None):
            commands.append((stage, command, sample))

    monkeypatch.setattr(cli, "WorkflowRuntime", FakeRuntime)
    result = CliRunner().invoke(
        cli.cli,
        [
            "sync-genome-views",
            "--db",
            str(tmp_path / "metatrawl.duckdb"),
            "--compare-dir",
            str(compare_dir),
            "--view-dir",
            str(tmp_path / "views"),
            "--workflow-config",
            str(workflow_config),
            "--neighbor-k",
            "12",
        ],
    )

    assert result.exit_code == 0, result.output
    assert [stage for stage, _, _ in commands] == ["genome_view", "genome_view"]
    assert sorted(sample for _, _, sample in commands) == ["genome_a", "genome_b"]
    for _, command, _ in commands:
        assert command[:2] == ["metatrawl", "sync-genome-views"]
        assert command[command.index("--neighbor-k") + 1] == "12"
        assert command[command.index("--max-null-fraction") + 1] == "0.35"
        assert "--max-null-samples" not in command
        assert "--local-child" in command
    assert "ready=2" in result.output

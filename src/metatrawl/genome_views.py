"""Build self-contained, per-genome visualization bundles."""

from __future__ import annotations

from dataclasses import dataclass
from datetime import date, datetime, timezone
import gzip
import json
import math
import os
from pathlib import Path
import re
from tempfile import TemporaryDirectory
from typing import Any
import uuid

import duckdb
import h5py
import numpy as np
import polars as pl

from metatrawl.logging import WorkflowLogger


VIEW_SCHEMA_VERSION = 2
VIEW_ARTIFACTS = (
    "samples.json",
    "sample_stats.json",
    "clusters.json",
    "similarity_ani.condensed.f32.gz",
    "total_positions.condensed.u64.gz",
    "clustermap.png",
    "dendrogram.svg",
    "dendrogram.json",
    "neighbor_network.json",
    "sample_stats.parquet",
    "distributions.json",
    "view_data.h5",
    "manifest.json",
)
SUPPORTED_LINKAGE_METHODS = ("single", "complete", "average", "weighted")


@dataclass(frozen=True)
class GenomeViewOptions:
    """Parameters controlling comparison filtering, clustering, and rendering."""

    min_comp_len: int = 10_000
    impute_ani: float = 97.0
    max_null_samples: int = 500
    linkage_method: str = "average"
    neighbor_k: int = 20
    clonal_cluster_threshold: float = 99.93
    strain_cluster_threshold: float = 99.8

    def validate(self) -> None:
        if self.min_comp_len < 0:
            raise ValueError("min_comp_len must be non-negative")
        if not 0 <= self.impute_ani <= 100:
            raise ValueError("impute_ani must be between 0 and 100")
        if self.max_null_samples < 0:
            raise ValueError("max_null_samples must be non-negative")
        if self.linkage_method not in SUPPORTED_LINKAGE_METHODS:
            choices = ", ".join(SUPPORTED_LINKAGE_METHODS)
            raise ValueError(f"linkage_method must be one of: {choices}")
        if self.neighbor_k < 1:
            raise ValueError("neighbor_k must be positive")
        for name, value in (
            ("clonal_cluster_threshold", self.clonal_cluster_threshold),
            ("strain_cluster_threshold", self.strain_cluster_threshold),
        ):
            if not 0 <= value <= 100:
                raise ValueError(f"{name} must be between 0 and 100")

    def as_dict(self) -> dict[str, object]:
        return {
            "min_comp_len": self.min_comp_len,
            "impute_ani": self.impute_ani,
            "max_null_samples": self.max_null_samples,
            "linkage_method": self.linkage_method,
            "neighbor_k": self.neighbor_k,
            "clonal_cluster_threshold": self.clonal_cluster_threshold,
            "strain_cluster_threshold": self.strain_cluster_threshold,
        }


@dataclass(frozen=True)
class GenomeViewSyncSummary:
    """Result from converging a collection of per-genome view bundles."""

    genomes: int
    ready: int
    generated: int
    up_to_date: int
    skipped: int
    failed: int


@dataclass(frozen=True)
class _PreparedGenomeView:
    genome: str
    samples: list[str]
    similarity_matrix: np.ndarray
    total_positions_matrix: np.ndarray
    null_fraction: np.ndarray
    linkage_matrix: np.ndarray
    leaf_order: np.ndarray
    clonal_clusters: np.ndarray
    strain_clusters: np.ndarray
    neighbor_edges: list[dict[str, object]]
    pair_ani: np.ndarray
    pair_total_positions: np.ndarray


@dataclass(frozen=True)
class _ComparisonDatabaseSchema:
    ani_column: str
    has_sample_catalog: bool
    has_genome_catalog: bool
    has_completed_pairs: bool

    @property
    def is_legacy(self) -> bool:
        return (
            self.ani_column != "genome_ani"
            or not self.has_sample_catalog
            or not self.has_genome_catalog
            or not self.has_completed_pairs
        )


class IncompleteComparisonError(RuntimeError):
    """Raised when a compare database is not ready for a stable view snapshot."""


def discover_compare_databases(compare_dir: Path, genomes: list[str] | None = None) -> list[Path]:
    """Return per-genome comparison DuckDB files, optionally restricted by name."""
    compare_dir = Path(compare_dir)
    if not compare_dir.is_dir():
        raise FileNotFoundError(f"Comparison directory does not exist: {compare_dir}")
    if genomes:
        paths = [compare_dir / f"{_safe_file_stem(genome)}.duckdb" for genome in _dedupe(genomes)]
        missing = [str(path) for path in paths if not path.is_file()]
        if missing:
            raise FileNotFoundError(f"Comparison database does not exist: {missing[0]}")
        return paths
    return sorted(path for path in compare_dir.iterdir() if path.is_file() and path.suffix.lower() == ".duckdb")


def sync_genome_views(
    *,
    db_file: Path,
    compare_dir: Path,
    view_dir: Path,
    genomes: list[str] | None = None,
    options: GenomeViewOptions | None = None,
    force: bool = False,
    logger: WorkflowLogger | None = None,
) -> GenomeViewSyncSummary:
    """Create or refresh one self-contained view bundle per comparison database."""
    logger = logger or WorkflowLogger()
    options = options or GenomeViewOptions()
    options.validate()
    compare_files = discover_compare_databases(compare_dir, genomes)
    view_dir = Path(view_dir)
    view_dir.mkdir(parents=True, exist_ok=True)
    logger.emit(step="sync-genome-views", status="start", genomes=len(compare_files), view_dir=view_dir)

    generated = up_to_date = skipped = failed = 0
    with duckdb.connect(str(Path(db_file)), read_only=True) as project_conn:
        for compare_file in compare_files:
            label = compare_file.stem
            logger.emit(step="sync-genome-views", status="genome-start", genome=label, compare_db=compare_file)
            try:
                status, genome = build_genome_view(
                    project_conn=project_conn,
                    compare_db=compare_file,
                    view_dir=view_dir,
                    options=options,
                    force=force,
                    logger=logger,
                )
            except IncompleteComparisonError as exc:
                skipped += 1
                logger.emit(step="sync-genome-views", status="skipped-incomplete", genome=label, error=exc)
            except duckdb.Error as exc:
                if "lock" in str(exc).lower():
                    skipped += 1
                    logger.emit(step="sync-genome-views", status="skipped-locked", genome=label, error=exc)
                    continue
                failed += 1
                logger.emit(step="sync-genome-views", status="failed", genome=label, error=exc)
            except (OSError, RuntimeError, ValueError) as exc:
                failed += 1
                logger.emit(step="sync-genome-views", status="failed", genome=label, error=exc)
            else:
                if status == "up-to-date":
                    up_to_date += 1
                else:
                    generated += 1
                logger.emit(step="sync-genome-views", status=status, genome=genome, output=view_dir / _safe_file_stem(genome))

    ready = generated + up_to_date
    write_view_catalog(view_dir)
    logger.emit(
        step="sync-genome-views",
        status="done",
        genomes=len(compare_files),
        ready=ready,
        generated=generated,
        up_to_date=up_to_date,
        skipped=skipped,
        failed=failed,
    )
    return GenomeViewSyncSummary(
        genomes=len(compare_files),
        ready=ready,
        generated=generated,
        up_to_date=up_to_date,
        skipped=skipped,
        failed=failed,
    )


def write_view_catalog(view_dir: Path) -> Path:
    """Atomically index all complete genome bundles for database-free discovery."""
    view_dir = Path(view_dir)
    view_dir.mkdir(parents=True, exist_ok=True)
    genomes = []
    for manifest_path in sorted(view_dir.glob("*/manifest.json")):
        try:
            manifest = json.loads(manifest_path.read_text())
        except (OSError, json.JSONDecodeError):
            continue
        if manifest.get("schema_version") != VIEW_SCHEMA_VERSION:
            continue
        genomes.append(
            {
                "genome": manifest.get("genome"),
                "path": manifest_path.parent.name,
                "manifest": f"{manifest_path.parent.name}/manifest.json",
                "sample_count": manifest.get("sample_count"),
                "neighbor_edge_count": manifest.get("neighbor_edge_count"),
                "clonal_cluster_count": manifest.get("clonal_cluster_count"),
                "strain_cluster_count": manifest.get("strain_cluster_count"),
                "generated_at": manifest.get("generated_at"),
                "source_signature": manifest.get("source_signature"),
            }
        )
    payload = {
        "schema_version": VIEW_SCHEMA_VERSION,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "genome_count": len(genomes),
        "genomes": genomes,
    }
    output = view_dir / "catalog.json"
    temporary = view_dir / f".catalog-{os.getpid()}-{uuid.uuid4().hex}.tmp"
    try:
        _write_json(temporary, payload)
        os.replace(temporary, output)
    finally:
        temporary.unlink(missing_ok=True)
    return output


def build_genome_view(
    *,
    project_conn: duckdb.DuckDBPyConnection,
    compare_db: Path,
    view_dir: Path,
    options: GenomeViewOptions,
    force: bool = False,
    logger: WorkflowLogger | None = None,
) -> tuple[str, str]:
    """Build one genome bundle and return ``(status, genome)``."""
    logger = logger or WorkflowLogger()
    options.validate()
    compare_db = Path(compare_db).expanduser().resolve()
    if not compare_db.is_file():
        raise FileNotFoundError(f"Comparison database does not exist: {compare_db}")

    with duckdb.connect(str(compare_db), read_only=True) as compare_conn:
        compare_schema = _detect_comparison_schema(compare_conn)
        genome = _single_genome(compare_conn, schema=compare_schema)
        sample_count, completed_pairs, expected_pairs, completion_source = _comparison_progress(
            compare_conn,
            schema=compare_schema,
        )
        if compare_schema.is_legacy:
            logger.emit(
                step="sync-genome-views",
                status="legacy-compare-schema",
                genome=genome,
                ani_column=compare_schema.ani_column,
                completion_source=completion_source,
            )
        if sample_count < 2:
            raise IncompleteComparisonError(
                f"genome={genome} needs at least two compared samples; found {sample_count}"
            )
        if completed_pairs < expected_pairs:
            raise IncompleteComparisonError(
                f"genome={genome} completed_pairs={completed_pairs} expected_pairs={expected_pairs}"
            )

        signature = _source_signature(
            compare_db,
            genome=genome,
            sample_count=sample_count,
            completed_pairs=completed_pairs,
            options=options,
            project_stats=_project_stats_signature(project_conn, genome=genome),
            comparison_schema={
                "ani_column": compare_schema.ani_column,
                "sample_catalog": compare_schema.has_sample_catalog,
                "genome_catalog": compare_schema.has_genome_catalog,
                "completed_pairs_table": compare_schema.has_completed_pairs,
                "completion_source": completion_source,
            },
        )
        output_dir = Path(view_dir) / _safe_file_stem(genome)
        if not force and _bundle_is_current(output_dir, signature):
            return "up-to-date", genome

        logger.emit(step="sync-genome-views", status="loading-comparisons", genome=genome)
        prepared = _prepare_genome_view(
            compare_conn,
            genome=genome,
            options=options,
            ani_column=compare_schema.ani_column,
        )

    logger.emit(step="sync-genome-views", status="loading-stats", genome=genome, samples=len(prepared.samples))
    sample_stats = _sample_stats(project_conn, genome=genome, samples=prepared.samples, leaf_order=prepared.leaf_order)
    _write_bundle(
        output_dir=output_dir,
        prepared=prepared,
        sample_stats=sample_stats,
        signature=signature,
        options=options,
    )
    return "generated", genome


def _prepare_genome_view(
    conn: duckdb.DuckDBPyConnection,
    *,
    genome: str,
    options: GenomeViewOptions,
    ani_column: str = "genome_ani",
) -> _PreparedGenomeView:
    from scipy.cluster.hierarchy import fcluster, leaves_list, linkage
    from scipy.spatial.distance import squareform

    quoted_ani = _quote_identifier(ani_column)
    pairs = conn.execute(
        f"""
        SELECT sample_1, sample_2,
               avg({quoted_ani}) AS genome_ani,
               max(total_positions) AS total_positions
        FROM matrix_compare_results
        WHERE genome = ?
          AND total_positions > ?
          AND {quoted_ani} IS NOT NULL
        GROUP BY sample_1, sample_2
        """,
        [genome, options.min_comp_len],
    ).pl()
    if pairs.is_empty():
        raise ValueError(
            f"No comparisons for genome={genome} pass min_comp_len={options.min_comp_len}"
        )

    samples = sorted(
        set(pairs.get_column("sample_1").to_list())
        | set(pairs.get_column("sample_2").to_list())
    )
    comparable_counts: dict[str, int] = {sample: 1 for sample in samples}
    for sample in pairs.get_column("sample_1"):
        comparable_counts[str(sample)] += 1
    for sample in pairs.get_column("sample_2"):
        comparable_counts[str(sample)] += 1
    total_sample_count = len(samples)
    excluded = {
        sample
        for sample, count in comparable_counts.items()
        if total_sample_count - count > options.max_null_samples
    }
    if excluded:
        samples = [sample for sample in samples if sample not in excluded]
        pairs = pairs.filter(
            (~pl.col("sample_1").is_in(excluded))
            & (~pl.col("sample_2").is_in(excluded))
        )
    if len(samples) < 2:
        raise ValueError(
            "At least two sufficiently connected samples are required. "
            "Relax max_null_samples or min_comp_len."
        )

    sample_index = {sample: idx for idx, sample in enumerate(samples)}
    similarity = np.full((len(samples), len(samples)), np.nan, dtype=np.float32)
    total_positions_matrix = np.zeros((len(samples), len(samples)), dtype=np.uint64)
    np.fill_diagonal(similarity, 100.0)
    idx_1 = np.fromiter(
        (sample_index[str(value)] for value in pairs.get_column("sample_1")),
        dtype=np.int64,
        count=pairs.height,
    )
    idx_2 = np.fromiter(
        (sample_index[str(value)] for value in pairs.get_column("sample_2")),
        dtype=np.int64,
        count=pairs.height,
    )
    ani_values = pairs.get_column("genome_ani").cast(pl.Float64).to_numpy()
    total_position_values = pairs.get_column("total_positions").cast(pl.UInt64).to_numpy()
    similarity[idx_1, idx_2] = ani_values
    similarity[idx_2, idx_1] = ani_values
    total_positions_matrix[idx_1, idx_2] = total_position_values
    total_positions_matrix[idx_2, idx_1] = total_position_values

    null_fraction = np.isnan(similarity).sum(axis=1).astype(np.float32) / len(samples)
    neighbor_edges = _top_neighbor_edges(
        similarity,
        total_positions_matrix=total_positions_matrix,
        samples=samples,
        k=options.neighbor_k,
    )
    np.nan_to_num(similarity, copy=False, nan=np.float32(options.impute_ani))

    distance_matrix = np.empty_like(similarity)
    np.divide(similarity, np.float32(100.0), out=distance_matrix)
    np.subtract(np.float32(1.0), distance_matrix, out=distance_matrix)
    np.fill_diagonal(distance_matrix, 0.0)
    linkage_matrix = linkage(squareform(distance_matrix, checks=False), method=options.linkage_method)
    leaf_order = leaves_list(linkage_matrix).astype(np.int32, copy=False)
    clonal_clusters = fcluster(
        linkage_matrix,
        t=1 - options.clonal_cluster_threshold / 100,
        criterion="distance",
    ).astype(np.int32, copy=False)
    strain_clusters = fcluster(
        linkage_matrix,
        t=1 - options.strain_cluster_threshold / 100,
        criterion="distance",
    ).astype(np.int32, copy=False)

    return _PreparedGenomeView(
        genome=genome,
        samples=samples,
        similarity_matrix=similarity,
        total_positions_matrix=total_positions_matrix,
        null_fraction=null_fraction,
        linkage_matrix=linkage_matrix,
        leaf_order=leaf_order,
        clonal_clusters=clonal_clusters,
        strain_clusters=strain_clusters,
        neighbor_edges=neighbor_edges,
        pair_ani=ani_values.astype(np.float32, copy=False),
        pair_total_positions=total_position_values,
    )


def _top_neighbor_edges(
    similarity: np.ndarray,
    *,
    total_positions_matrix: np.ndarray,
    samples: list[str],
    k: int,
) -> list[dict[str, object]]:
    selected: dict[tuple[int, int], float] = {}
    for source_idx in range(len(samples)):
        row = similarity[source_idx]
        candidates = np.flatnonzero(np.isfinite(row))
        candidates = candidates[candidates != source_idx]
        if candidates.size == 0:
            continue
        values = row[candidates]
        keep = min(k, candidates.size)
        if keep < candidates.size:
            positions = np.argpartition(-values, keep - 1)[:keep]
            candidates = candidates[positions]
            values = values[positions]
        order = np.argsort(-values, kind="stable")
        for target_idx, ani in zip(candidates[order], values[order], strict=True):
            edge = (min(source_idx, int(target_idx)), max(source_idx, int(target_idx)))
            selected[edge] = max(selected.get(edge, float("-inf")), float(ani))
    return [
        {
            "source": samples[source],
            "target": samples[target],
            "source_index": source,
            "target_index": target,
            "ani": ani,
            "total_positions": int(total_positions_matrix[source, target]),
        }
        for (source, target), ani in sorted(selected.items(), key=lambda item: (-item[1], item[0]))
    ]


def _sample_stats(
    conn: duckdb.DuckDBPyConnection,
    *,
    genome: str,
    samples: list[str],
    leaf_order: np.ndarray,
) -> pl.DataFrame:
    stats = conn.execute(
        """
        WITH abundance AS (
          SELECT sample_id, max(abundance) AS sylph_abundance
          FROM sylph_abundance
          WHERE genome = ? OR accession = ?
          GROUP BY sample_id
        )
        SELECT gs.* EXCLUDE (genome), a.sylph_abundance
        FROM genome_stats gs
        LEFT JOIN abundance a ON a.sample_id = gs.sample_id
        WHERE gs.genome = ?
        """,
        [genome, genome, genome],
    ).pl()
    leaf_rank = np.empty(len(samples), dtype=np.int32)
    leaf_rank[leaf_order] = np.arange(len(samples), dtype=np.int32)
    result = (
        pl.DataFrame({"sample_id": samples, "leaf_order": leaf_rank})
        .join(stats, on="sample_id", how="left")
        .sort("leaf_order")
    )
    if result.height != len(samples):
        raise ValueError(
            f"Expected at most one genome_stats row per sample for genome={genome}; "
            f"found {result.height} rows for {len(samples)} samples"
        )
    return result


def _write_bundle(
    *,
    output_dir: Path,
    prepared: _PreparedGenomeView,
    sample_stats: pl.DataFrame,
    signature: dict[str, object],
    options: GenomeViewOptions,
) -> None:
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = output_dir / "manifest.json"
    manifest_path.unlink(missing_ok=True)
    with TemporaryDirectory(prefix=f".{output_dir.name}-", dir=output_dir.parent) as temp_name:
        staging = Path(temp_name)
        sample_stats.write_parquet(staging / "sample_stats.parquet")
        _write_samples_json(staging / "samples.json", prepared=prepared)
        _write_sample_stats_json(
            staging / "sample_stats.json",
            genome=prepared.genome,
            sample_stats=sample_stats,
        )
        _write_clusters_json(staging / "clusters.json", prepared=prepared, options=options)
        _write_condensed_matrix_gzip(
            staging / "similarity_ani.condensed.f32.gz",
            prepared.similarity_matrix,
            dtype="<f4",
        )
        _write_condensed_matrix_gzip(
            staging / "total_positions.condensed.u64.gz",
            prepared.total_positions_matrix,
            dtype="<u8",
        )
        _write_view_hdf5(staging / "view_data.h5", prepared=prepared, options=options)
        _write_dendrogram_json(staging / "dendrogram.json", prepared=prepared, options=options)
        _write_network_json(
            staging / "neighbor_network.json",
            prepared=prepared,
            sample_stats=sample_stats,
            options=options,
        )
        _write_distributions_json(
            staging / "distributions.json",
            prepared=prepared,
            sample_stats=sample_stats,
        )
        _render_clustermap(staging / "clustermap.png", prepared=prepared)
        _render_dendrogram(staging / "dendrogram.svg", prepared=prepared, options=options)

        manifest = {
            "schema_version": VIEW_SCHEMA_VERSION,
            "genome": prepared.genome,
            "generated_at": datetime.now(timezone.utc).isoformat(),
            "sample_count": len(prepared.samples),
            "neighbor_edge_count": len(prepared.neighbor_edges),
            "clonal_cluster_count": int(np.unique(prepared.clonal_clusters).size),
            "strain_cluster_count": int(np.unique(prepared.strain_clusters).size),
            "source_signature": signature,
            "options": options.as_dict(),
            "artifacts": list(VIEW_ARTIFACTS),
            "files": _artifact_descriptors(staging, prepared=prepared),
        }
        _write_json(staging / "manifest.json", manifest)

        for artifact in VIEW_ARTIFACTS:
            if artifact == "manifest.json":
                continue
            os.replace(staging / artifact, output_dir / artifact)
        os.replace(staging / "manifest.json", manifest_path)


def _write_samples_json(path: Path, *, prepared: _PreparedGenomeView) -> None:
    leaf_rank = np.empty(len(prepared.samples), dtype=np.int32)
    leaf_rank[prepared.leaf_order] = np.arange(len(prepared.samples), dtype=np.int32)
    records = [
        {
            "sample_index": sample_index,
            "sample_id": sample,
            "leaf_order": int(leaf_rank[sample_index]),
            "null_fraction": float(prepared.null_fraction[sample_index]),
            "clonal_cluster": int(prepared.clonal_clusters[sample_index]),
            "strain_cluster": int(prepared.strain_clusters[sample_index]),
        }
        for sample_index, sample in enumerate(prepared.samples)
    ]
    _write_json(
        path,
        {
            "schema_version": VIEW_SCHEMA_VERSION,
            "genome": prepared.genome,
            "sample_count": len(prepared.samples),
            "sample_order": "matrix-axis",
            "samples": records,
            "leaf_order": prepared.leaf_order.astype(int).tolist(),
            "ordered_samples": [prepared.samples[index] for index in prepared.leaf_order],
        },
    )


def _write_sample_stats_json(path: Path, *, genome: str, sample_stats: pl.DataFrame) -> None:
    columns = {
        column: [_json_scalar(value) for value in sample_stats.get_column(column).to_list()]
        for column in sample_stats.columns
    }
    schema = [
        {
            "name": column,
            "dtype": str(sample_stats.schema[column]),
            "nullable": sample_stats.get_column(column).null_count() > 0,
        }
        for column in sample_stats.columns
    ]
    _write_json(
        path,
        {
            "schema_version": VIEW_SCHEMA_VERSION,
            "genome": genome,
            "orientation": "columnar",
            "row_count": sample_stats.height,
            "schema": schema,
            "columns": columns,
        },
    )


def _write_clusters_json(
    path: Path,
    *,
    prepared: _PreparedGenomeView,
    options: GenomeViewOptions,
) -> None:
    assignments = [
        {
            "sample_index": sample_index,
            "sample_id": sample,
            "clonal_cluster": int(prepared.clonal_clusters[sample_index]),
            "strain_cluster": int(prepared.strain_clusters[sample_index]),
        }
        for sample_index, sample in enumerate(prepared.samples)
    ]
    _write_json(
        path,
        {
            "schema_version": VIEW_SCHEMA_VERSION,
            "genome": prepared.genome,
            "method": {
                "distance": "1 - genome_ani / 100",
                "linkage": options.linkage_method,
            },
            "thresholds": {
                "clonal_ani": options.clonal_cluster_threshold,
                "strain_ani": options.strain_cluster_threshold,
            },
            "assignments": assignments,
            "clusters": {
                "clonal": _cluster_groups(
                    prepared.samples,
                    prepared.clonal_clusters,
                    leaf_order=prepared.leaf_order,
                ),
                "strain": _cluster_groups(
                    prepared.samples,
                    prepared.strain_clusters,
                    leaf_order=prepared.leaf_order,
                ),
            },
        },
    )


def _cluster_groups(
    samples: list[str],
    labels: np.ndarray,
    *,
    leaf_order: np.ndarray,
) -> list[dict[str, object]]:
    groups: dict[int, list[dict[str, object]]] = {}
    for sample_index in leaf_order:
        index = int(sample_index)
        cluster = int(labels[index])
        groups.setdefault(cluster, []).append(
            {"sample_index": index, "sample_id": samples[index]}
        )
    return [
        {
            "cluster_id": cluster,
            "sample_count": len(members),
            "members": members,
        }
        for cluster, members in sorted(groups.items())
    ]


def _write_condensed_matrix_gzip(path: Path, matrix: np.ndarray, *, dtype: str) -> None:
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError("Condensed matrix export requires a square matrix")
    output_dtype = np.dtype(dtype)
    with path.open("wb") as raw_handle:
        with gzip.GzipFile(
            filename="",
            mode="wb",
            fileobj=raw_handle,
            compresslevel=6,
            mtime=0,
        ) as gzip_handle:
            for row_index in range(matrix.shape[0] - 1):
                row = np.ascontiguousarray(
                    matrix[row_index, row_index + 1 :],
                    dtype=output_dtype,
                )
                gzip_handle.write(memoryview(row).cast("B"))


def _artifact_descriptors(
    staging: Path,
    *,
    prepared: _PreparedGenomeView,
) -> dict[str, dict[str, object]]:
    sample_count = len(prepared.samples)
    definitions: dict[str, dict[str, object]] = {
        "samples": {
            "path": "samples.json",
            "format": "json",
            "media_type": "application/json",
            "description": "Sample indices, leaf order, null fraction, and cluster assignments.",
        },
        "sample_stats": {
            "path": "sample_stats.json",
            "format": "columnar-json",
            "media_type": "application/json",
            "description": "Browser-ready sample statistics in dendrogram leaf order.",
        },
        "sample_stats_parquet": {
            "path": "sample_stats.parquet",
            "format": "parquet",
            "media_type": "application/vnd.apache.parquet",
            "description": "Typed sample statistics for Python, R, DuckDB, and Polars.",
        },
        "clusters": {
            "path": "clusters.json",
            "format": "json",
            "media_type": "application/json",
            "description": "Reusable clonal and strain cluster assignments and memberships.",
        },
        "dendrogram": {
            "path": "dendrogram.json",
            "format": "scipy-linkage-json",
            "media_type": "application/json",
            "description": "Linkage matrix and explicit merge-tree topology.",
        },
        "neighbor_network": {
            "path": "neighbor_network.json",
            "format": "node-link-json",
            "media_type": "application/json",
            "description": "Sparse top-neighbor sample graph.",
        },
        "distributions": {
            "path": "distributions.json",
            "format": "histogram-json",
            "media_type": "application/json",
            "description": "Precomputed ANI, overlap, and sample-statistic histograms.",
        },
        "similarity_matrix": {
            "path": "similarity_ani.condensed.f32.gz",
            "format": "scipy-condensed-upper-triangle",
            "media_type": "application/gzip",
            "compression": "gzip",
            "dtype": "float32",
            "byte_order": "little",
            "order": "C",
            "shape": [sample_count * (sample_count - 1) // 2],
            "matrix_shape": [sample_count, sample_count],
            "axis_order": "samples.json samples[].sample_index",
            "indexing": "k = n*i - i*(i+1)/2 + j-i-1 for i < j",
            "diagonal_value": 100.0,
            "description": "Symmetric ANI percentage matrix with missing cells imputed.",
        },
        "total_positions_matrix": {
            "path": "total_positions.condensed.u64.gz",
            "format": "scipy-condensed-upper-triangle",
            "media_type": "application/gzip",
            "compression": "gzip",
            "dtype": "uint64",
            "byte_order": "little",
            "order": "C",
            "shape": [sample_count * (sample_count - 1) // 2],
            "matrix_shape": [sample_count, sample_count],
            "axis_order": "samples.json samples[].sample_index",
            "indexing": "k = n*i - i*(i+1)/2 + j-i-1 for i < j",
            "diagonal_value": 0,
            "missing_value": 0,
            "description": "Pairwise compared-position counts; zero identifies imputed ANI cells and the diagonal.",
        },
        "analysis_hdf5": {
            "path": "view_data.h5",
            "format": "hdf5",
            "media_type": "application/x-hdf5",
            "description": "Reusable matrices, linkage, ordering, and assignments.",
        },
        "clustermap_preview": {
            "path": "clustermap.png",
            "format": "png",
            "media_type": "image/png",
            "description": "Static preview; the web heatmap should use the binary matrices.",
        },
        "dendrogram_preview": {
            "path": "dendrogram.svg",
            "format": "svg",
            "media_type": "image/svg+xml",
            "description": "Static preview; the web dendrogram should use dendrogram.json.",
        },
    }
    for descriptor in definitions.values():
        descriptor["size_bytes"] = (staging / str(descriptor["path"])).stat().st_size
    return definitions


def _write_view_hdf5(path: Path, *, prepared: _PreparedGenomeView, options: GenomeViewOptions) -> None:
    string_dtype = h5py.string_dtype(encoding="utf-8")
    sample_count = len(prepared.samples)
    chunk_size = max(1, min(sample_count, 512))
    with h5py.File(path, "w") as handle:
        handle.attrs["schema_version"] = VIEW_SCHEMA_VERSION
        handle.attrs["genome"] = prepared.genome
        for key, value in options.as_dict().items():
            handle.attrs[key] = value
        handle.create_dataset("samples", data=np.asarray(prepared.samples, dtype=object), dtype=string_dtype)
        handle.create_dataset(
            "similarity_ani",
            data=prepared.similarity_matrix,
            compression="gzip",
            compression_opts=4,
            shuffle=True,
            chunks=(chunk_size, chunk_size),
        )
        handle.create_dataset(
            "total_positions",
            data=prepared.total_positions_matrix,
            compression="gzip",
            compression_opts=4,
            shuffle=True,
            chunks=(chunk_size, chunk_size),
        )
        handle.create_dataset("linkage", data=prepared.linkage_matrix, compression="gzip", shuffle=True)
        handle.create_dataset("leaf_order", data=prepared.leaf_order)
        handle.create_dataset("null_fraction", data=prepared.null_fraction)
        handle.create_dataset("clonal_cluster", data=prepared.clonal_clusters)
        handle.create_dataset("strain_cluster", data=prepared.strain_clusters)


def _write_dendrogram_json(path: Path, *, prepared: _PreparedGenomeView, options: GenomeViewOptions) -> None:
    sample_count = len(prepared.samples)
    merges = [
        {
            "id": sample_count + merge_index,
            "left": int(row[0]),
            "right": int(row[1]),
            "distance": float(row[2]),
            "size": int(row[3]),
        }
        for merge_index, row in enumerate(prepared.linkage_matrix)
    ]
    _write_json(
        path,
        {
            "genome": prepared.genome,
            "distance": "1 - genome_ani / 100",
            "linkage_method": options.linkage_method,
            "samples": prepared.samples,
            "leaf_order": prepared.leaf_order.astype(int).tolist(),
            "ordered_samples": [prepared.samples[idx] for idx in prepared.leaf_order],
            "linkage": {
                "format": "scipy",
                "columns": ["left", "right", "distance", "size"],
                "rows": prepared.linkage_matrix.tolist(),
            },
            "tree": {
                "leaf_ids": list(range(sample_count)),
                "merges": merges,
                "root_id": sample_count * 2 - 2,
            },
            "clonal_cluster_threshold": options.clonal_cluster_threshold,
            "strain_cluster_threshold": options.strain_cluster_threshold,
        },
    )


def _write_network_json(
    path: Path,
    *,
    prepared: _PreparedGenomeView,
    sample_stats: pl.DataFrame,
    options: GenomeViewOptions,
) -> None:
    stats_by_sample = {
        str(row["sample_id"]): row
        for row in sample_stats.to_dicts()
    }
    nodes = []
    for idx, sample in enumerate(prepared.samples):
        stats = stats_by_sample.get(sample, {})
        nodes.append(
            {
                "id": sample,
                "leaf_order": _json_scalar(stats.get("leaf_order")),
                "clonal_cluster": int(prepared.clonal_clusters[idx]),
                "strain_cluster": int(prepared.strain_clusters[idx]),
                "coverage": _json_scalar(stats.get("coverage")),
                "breadth": _json_scalar(stats.get("breadth")),
                "ber": _json_scalar(stats.get("ber")),
                "ref_ani": _json_scalar(stats.get("ref_ani")),
                "sylph_abundance": _json_scalar(stats.get("sylph_abundance")),
            }
        )
    _write_json(
        path,
        {
            "genome": prepared.genome,
            "neighbor_k": options.neighbor_k,
            "selection": "union of each sample's top-k comparable ANI neighbors",
            "nodes": nodes,
            "edges": prepared.neighbor_edges,
        },
    )


def _write_distributions_json(
    path: Path,
    *,
    prepared: _PreparedGenomeView,
    sample_stats: pl.DataFrame,
) -> None:
    histograms: dict[str, object] = {
        "genome_ani": _histogram(prepared.pair_ani, bins=250, value_range=(95.0, 100.0)),
        "genome_ani_full": _histogram(prepared.pair_ani, bins=200, value_range=(0.0, 100.0)),
        "total_positions": _histogram(prepared.pair_total_positions, bins=100),
    }
    for column in ("coverage", "breadth", "ber", "ref_ani", "sylph_abundance"):
        if column not in sample_stats.columns:
            continue
        values = sample_stats.get_column(column).drop_nulls().cast(pl.Float64).to_numpy()
        histograms[column] = _histogram(values, bins=100)
    _write_json(path, {"genome": prepared.genome, "histograms": histograms})


def _render_clustermap(path: Path, *, prepared: _PreparedGenomeView) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import seaborn as sns

    sample_count = len(prepared.samples)
    figure_size = min(30.0, max(10.0, 8.0 + math.log10(max(sample_count, 2)) * 3.0))
    show_labels = sample_count <= 100
    labels: list[str] | bool = prepared.samples if show_labels else False
    grid = sns.clustermap(
        prepared.similarity_matrix,
        row_linkage=prepared.linkage_matrix,
        col_linkage=prepared.linkage_matrix,
        cmap="rocket",
        vmin=max(0.0, float(np.nanpercentile(prepared.similarity_matrix, 1))),
        vmax=100.0,
        xticklabels=labels,
        yticklabels=labels,
        figsize=(figure_size, figure_size),
        cbar_kws={"label": "PopANI (%)"},
    )
    grid.fig.suptitle(f"{prepared.genome} clustermap", y=1.01)
    grid.ax_heatmap.set_xlabel("Samples")
    grid.ax_heatmap.set_ylabel("Samples")
    grid.fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(grid.fig)


def _render_dendrogram(path: Path, *, prepared: _PreparedGenomeView, options: GenomeViewOptions) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from scipy.cluster.hierarchy import dendrogram

    sample_count = len(prepared.samples)
    show_labels = sample_count <= 500
    figure_height = max(10.0, min(80.0, sample_count * 0.08 if show_labels else 14.0))
    fig, ax = plt.subplots(figsize=(12, figure_height))
    dendrogram(
        prepared.linkage_matrix,
        labels=prepared.samples if show_labels else None,
        orientation="left",
        no_labels=not show_labels,
        color_threshold=1 - options.strain_cluster_threshold / 100,
        above_threshold_color="#808080",
        leaf_font_size=7,
    )
    ax.axvline(1 - options.clonal_cluster_threshold / 100, color="#D1495B", linestyle="--", linewidth=1)
    ax.axvline(1 - options.strain_cluster_threshold / 100, color="#00798C", linestyle="--", linewidth=1)
    ax.set_title(f"{prepared.genome} dendrogram")
    ax.set_xlabel("PopANI distance (1 - ANI/100)")
    ax.set_ylabel("Samples")
    fig.tight_layout()
    fig.savefig(path, format="svg", bbox_inches="tight")
    plt.close(fig)


def _comparison_progress(
    conn: duckdb.DuckDBPyConnection,
    *,
    schema: _ComparisonDatabaseSchema,
) -> tuple[int, int, int, str]:
    if schema.has_sample_catalog:
        sample_count = int(conn.execute("SELECT count(*) FROM matrix_compare_samples").fetchone()[0])
    else:
        sample_count = int(
            conn.execute(
                """
                SELECT count(*)
                FROM (
                  SELECT sample_1 AS sample FROM matrix_compare_results
                  UNION
                  SELECT sample_2 AS sample FROM matrix_compare_results
                )
                """
            ).fetchone()[0]
        )
    if schema.has_genome_catalog:
        genome_count = int(conn.execute("SELECT count(*) FROM matrix_compare_genomes").fetchone()[0])
    else:
        genome_count = int(
            conn.execute(
                "SELECT count(DISTINCT genome) FROM matrix_compare_results"
            ).fetchone()[0]
        )
    if schema.has_completed_pairs:
        completed_pairs = int(
            conn.execute("SELECT count(*) FROM matrix_compare_completed_pair_genomes").fetchone()[0]
        )
        completion_source = "matrix_compare_completed_pair_genomes"
    else:
        completed_pairs = int(
            conn.execute(
                """
                SELECT count(*)
                FROM (
                  SELECT DISTINCT
                    CASE WHEN sample_1 <= sample_2 THEN sample_1 ELSE sample_2 END AS sample_1,
                    CASE WHEN sample_1 <= sample_2 THEN sample_2 ELSE sample_1 END AS sample_2,
                    genome
                  FROM matrix_compare_results
                )
                """
            ).fetchone()[0]
        )
        completion_source = "distinct_result_rows"
    expected_pairs = sample_count * (sample_count - 1) // 2 * genome_count
    return sample_count, completed_pairs, expected_pairs, completion_source


def _single_genome(
    conn: duckdb.DuckDBPyConnection,
    *,
    schema: _ComparisonDatabaseSchema,
) -> str:
    if schema.has_genome_catalog:
        rows = conn.execute(
            "SELECT genome FROM matrix_compare_genomes ORDER BY genome_idx"
        ).fetchall()
    else:
        rows = conn.execute(
            "SELECT DISTINCT genome FROM matrix_compare_results ORDER BY genome"
        ).fetchall()
    genomes = [str(row[0]) for row in rows]
    if len(genomes) != 1:
        raise ValueError(f"Genome views require one genome per comparison database; found {len(genomes)}")
    return genomes[0]


def _detect_comparison_schema(
    conn: duckdb.DuckDBPyConnection,
) -> _ComparisonDatabaseSchema:
    existing = {str(row[0]).lower() for row in conn.execute("SHOW TABLES").fetchall()}
    if "matrix_compare_results" not in existing:
        raise ValueError("Comparison database is missing required table: matrix_compare_results")
    result_columns = {
        str(row[0]).lower(): str(row[0])
        for row in conn.execute(
            """
            SELECT column_name
            FROM information_schema.columns
            WHERE table_schema = current_schema()
              AND table_name = 'matrix_compare_results'
            """
        ).fetchall()
    }
    required_columns = {"sample_1", "sample_2", "genome", "total_positions"}
    missing_columns = sorted(required_columns - set(result_columns))
    if missing_columns:
        raise ValueError(
            "Comparison result table is missing required column: "
            + missing_columns[0]
        )
    ani_column = next(
        (
            result_columns[candidate]
            for candidate in ("genome_ani", "genome_pop_ani")
            if candidate in result_columns
        ),
        None,
    )
    if ani_column is None:
        raise ValueError(
            "Comparison result table has no supported ANI column; "
            "expected genome_ani or genome_pop_ani"
        )
    return _ComparisonDatabaseSchema(
        ani_column=ani_column,
        has_sample_catalog="matrix_compare_samples" in existing,
        has_genome_catalog="matrix_compare_genomes" in existing,
        has_completed_pairs="matrix_compare_completed_pair_genomes" in existing,
    )


def _source_signature(
    compare_db: Path,
    *,
    genome: str,
    sample_count: int,
    completed_pairs: int,
    options: GenomeViewOptions,
    project_stats: dict[str, object],
    comparison_schema: dict[str, object],
) -> dict[str, object]:
    stat = compare_db.stat()
    return {
        "compare_db": str(compare_db),
        "size": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
        "genome": genome,
        "sample_count": sample_count,
        "completed_pairs": completed_pairs,
        "options": options.as_dict(),
        "project_stats": project_stats,
        "comparison_schema": comparison_schema,
    }


def _project_stats_signature(
    conn: duckdb.DuckDBPyConnection,
    *,
    genome: str,
) -> dict[str, object]:
    return {
        "genome_stats": _filtered_table_signature(
            conn,
            table="genome_stats",
            where_sql="genome = ?",
            parameters=[genome],
        ),
        "sylph_abundance": _filtered_table_signature(
            conn,
            table="sylph_abundance",
            where_sql="genome = ? OR accession = ?",
            parameters=[genome, genome],
        ),
    }


def _filtered_table_signature(
    conn: duckdb.DuckDBPyConnection,
    *,
    table: str,
    where_sql: str,
    parameters: list[object],
) -> dict[str, object]:
    columns = [
        str(row[0])
        for row in conn.execute(
            """
            SELECT column_name
            FROM information_schema.columns
            WHERE table_schema = current_schema()
              AND table_name = ?
            ORDER BY ordinal_position
            """,
            [table],
        ).fetchall()
    ]
    if not columns:
        return {"rows": 0, "checksum": "0"}
    hash_arguments = ", ".join(_quote_identifier(column) for column in columns)
    row_count, checksum = conn.execute(
        f"""
        SELECT count(*), COALESCE(bit_xor(hash({hash_arguments})), 0)
        FROM {_quote_identifier(table)}
        WHERE {where_sql}
        """,
        parameters,
    ).fetchone()
    return {"rows": int(row_count), "checksum": str(int(checksum))}


def _bundle_is_current(output_dir: Path, signature: dict[str, object]) -> bool:
    manifest = output_dir / "manifest.json"
    if not manifest.is_file() or any(not (output_dir / name).is_file() for name in VIEW_ARTIFACTS):
        return False
    try:
        payload = json.loads(manifest.read_text())
    except (OSError, json.JSONDecodeError):
        return False
    return payload.get("schema_version") == VIEW_SCHEMA_VERSION and payload.get("source_signature") == signature


def _histogram(
    values: np.ndarray,
    *,
    bins: int,
    value_range: tuple[float, float] | None = None,
) -> dict[str, object]:
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        return {"counts": [], "bin_edges": [], "value_count": 0}
    if value_range is None and float(finite.min()) == float(finite.max()):
        center = float(finite[0])
        padding = max(abs(center) * 0.01, 0.5)
        value_range = (center - padding, center + padding)
    counts, edges = np.histogram(finite, bins=bins, range=value_range)
    return {
        "counts": counts.astype(int).tolist(),
        "bin_edges": edges.astype(float).tolist(),
        "value_count": int(finite.size),
    }


def _write_json(path: Path, payload: dict[str, object]) -> None:
    path.write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n")


def _json_scalar(value: Any) -> object:
    if value is None:
        return None
    if isinstance(value, (float, np.floating)) and not math.isfinite(float(value)):
        return None
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, (date, datetime)):
        return value.isoformat()
    return value


def _safe_file_stem(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", value).strip("_") or "genome"


def _dedupe(values: list[str]) -> list[str]:
    return list(dict.fromkeys(value.strip() for value in values if value.strip()))


def _quote_identifier(value: str) -> str:
    return '"' + value.replace('"', '""') + '"'

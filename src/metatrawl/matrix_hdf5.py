"""Direct DuckDB-to-HDF5 matrix writer for MetaTrawl.

This module intentionally writes ZipStrain-compatible HDF5 matrix stores without
staging temporary profile parquet files. It mirrors the current ZipStrain HDF5
contract closely enough that existing ZipStrain matrix compare code can consume
files produced here.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import duckdb
import numpy as np
import polars as pl

from metatrawl import db

BuildProgressCallback = Callable[[dict[str, object]], None]

COUNT_DTYPES = {
    "uint16": np.uint16,
    "uint32": np.uint32,
}
MATRIX_BUILD_MIN_COV = 5
FILTERED_PRESENCE_MATRIX_VALUE_SEMANTICS = "allele_presence_after_cov_filter"
CURRENT_MATRIX_HDF5_LAYOUT = "per_genome_sample_major_dense_matrix_hdf5"
CURRENT_MATRIX_HDF5_SPARSE_LAYOUT = "per_genome_sample_major_sparse_indices_matrix_hdf5"
MATRIX_HDF5_FILE_VERSION = "1"
MATRIX_HDF5_CONTRACT_GENOMES_GROUP = "contract_genomes"
MATRIX_HDF5_CONTRACT_SCAFFOLDS_GROUP = "contract_genome_scaffolds"


@dataclass(frozen=True)
class ScaffoldSpec:
    scaffold_idx: int
    genome: str
    chrom: str
    index_base: int
    vector_length: int
    min_pos: int
    max_pos: int


@dataclass(frozen=True)
class GenomeSpec:
    genome_idx: int
    genome: str
    matrix_length: int
    true_length: int
    scaffold_count: int


@dataclass(frozen=True)
class GenomeScaffoldOffset:
    genome_idx: int
    scaffold_ordinal: int
    genome: str
    chrom: str
    axis_start: int
    axis_end: int
    index_base: int
    vector_length: int
    min_pos: int
    max_pos: int


@dataclass(frozen=True)
class GeneRangeSpec:
    gene_idx: int
    gene: str
    genome_idx: int
    genome: str
    chrom: str
    axis_start: int
    axis_end: int


@dataclass(frozen=True)
class DirectMatrixBuildSummary:
    output_file: Path
    sample_count: int
    stored_rows: int
    storage_layout: str
    genome_scope: str


@dataclass(frozen=True)
class DirectMatrixAppendSummary:
    output_file: Path
    appended_sample_count: int
    total_sample_count: int
    stored_rows: int


def build_matrix_hdf5_from_duckdb(
    conn: duckdb.DuckDBPyConnection,
    *,
    sample_ids: list[str],
    output_file: Path,
    genome: str,
    bed_file: Path,
    stb_file: Path,
    gene_range_table: Path | None = None,
    count_dtype: str = "uint16",
    memory_limit_gb: float = 16.0,
    export_batch_mb: float = 128.0,
    duckdb_threads: int = 1,
    sparse: bool = False,
    progress_callback: BuildProgressCallback | None = None,
) -> DirectMatrixBuildSummary:
    """Build a ZipStrain-compatible HDF5 matrix directly from MetaTrawl DuckDB rows."""
    if count_dtype not in COUNT_DTYPES:
        raise ValueError(f"Unsupported count dtype '{count_dtype}'. Choose one of {', '.join(COUNT_DTYPES)}.")
    if export_batch_mb <= 0:
        raise ValueError("export_batch_mb must be > 0")
    _configure_duckdb_for_matrix(conn, memory_limit_gb=memory_limit_gb, threads=duckdb_threads)

    genome_scope = None if genome == "all" else genome
    contract_scaffolds = _collect_scaffold_specs_from_bed_and_stb(
        bed_file=Path(bed_file),
        stb_file=Path(stb_file),
        genome=genome_scope,
    )
    contract_genomes, contract_genome_scaffolds = _build_genome_specs(contract_scaffolds)
    observed_genome_names = _observed_genomes_for_samples(conn, sample_ids=sample_ids, genome=genome_scope)
    if not observed_genome_names:
        detail = genome_scope if genome_scope is not None else "all"
        raise ValueError(f"No profile rows found for genome scope: {detail}")
    unknown_observed = observed_genome_names - {spec.genome for spec in contract_genomes}
    if unknown_observed:
        raise ValueError(
            "Profiles contain genomes that are missing from the provided BED/STB contract: "
            + ", ".join(sorted(unknown_observed))
        )
    genomes, genome_scaffolds = _subset_contract_to_genomes(
        contract_genomes=contract_genomes,
        contract_genome_scaffolds=contract_genome_scaffolds,
        genome_names=observed_genome_names,
    )
    gene_ranges = (
        _collect_gene_range_specs(gene_range_table=Path(gene_range_table), genome_scaffolds=genome_scaffolds)
        if gene_range_table is not None
        else []
    )

    output_file = Path(output_file).expanduser().resolve()
    output_file.parent.mkdir(parents=True, exist_ok=True)
    if output_file.exists():
        raise FileExistsError(f"Output file already exists: {output_file}")

    memory_limit_bytes = _memory_limit_bytes(memory_limit_gb)
    total_work = len(sample_ids) * len(genomes)
    completed_work = 0
    stored_rows = 0
    _emit_progress(progress_callback, "start", completed_work, total_work, "", genome_scope or "all", stored_rows)

    h5py = _import_h5py()
    tmp_output = output_file.with_suffix(output_file.suffix + ".tmp")
    if tmp_output.exists():
        tmp_output.unlink()
    build_succeeded = False
    try:
        sample_rows = [(idx, sample_id) for idx, sample_id in enumerate(sample_ids)]
        metadata_rows = _metadata_rows(
            genome=genome_scope or "all",
            count_dtype=count_dtype,
            memory_limit_gb=memory_limit_gb,
            export_batch_mb=export_batch_mb,
            sparse=sparse,
            bed_file=bed_file,
            stb_file=stb_file,
            gene_range_table=gene_range_table,
            has_gene_ranges=bool(gene_ranges),
        )
        scaffolds_by_genome_idx = _group_scaffolds_by_genome(genome_scaffolds)
        with h5py.File(str(tmp_output), "w") as handle:
            _initialize_hdf5_store(
                handle,
                h5py_module=h5py,
                metadata_rows=metadata_rows,
                sample_rows=sample_rows,
                genomes=genomes,
                genome_scaffolds=genome_scaffolds,
                contract_genomes=contract_genomes,
                contract_genome_scaffolds=contract_genome_scaffolds,
                gene_ranges=gene_ranges,
                count_dtype=count_dtype,
                export_batch_mb=export_batch_mb,
                sparse=sparse,
            )
            sparse_states = _initial_sparse_states(handle, genomes) if sparse else {}
            for sample_row, sample_id in sample_rows:
                for spec in genomes:
                    _check_matrix_memory(spec, count_dtype=count_dtype, memory_limit_bytes=memory_limit_bytes)
                    _emit_progress(progress_callback, "processing", completed_work, total_work, sample_id, spec.genome, stored_rows)
                    matrix = _load_duckdb_sample_genome_matrix(
                        conn,
                        sample_id=sample_id,
                        genome_spec=spec,
                        genome_offsets=scaffolds_by_genome_idx[spec.genome_idx],
                        count_dtype=count_dtype,
                    )
                    if sparse:
                        indptr_ds, indices_ds, current_nnz = sparse_states[spec.genome_idx]
                        current_nnz = _append_sparse_hdf5_matrix_row(
                            indptr_dataset=indptr_ds,
                            indices_dataset=indices_ds,
                            sample_row=sample_row,
                            flat_indices=_dense_matrix_to_sparse_flat_indices(matrix),
                            current_nnz=current_nnz,
                        )
                        sparse_states[spec.genome_idx] = (indptr_ds, indices_ds, current_nnz)
                    else:
                        handle[_hdf5_matrix_dataset_path(spec.genome_idx)][sample_row, :, :] = matrix
                    stored_rows += 1
                    completed_work += 1
                    _emit_progress(progress_callback, "advance", completed_work, total_work, sample_id, spec.genome, stored_rows)
        tmp_output.replace(output_file)
        build_succeeded = True
    finally:
        if not build_succeeded and tmp_output.exists():
            tmp_output.unlink(missing_ok=True)
        if not build_succeeded and output_file.exists():
            output_file.unlink(missing_ok=True)

    _emit_progress(progress_callback, "done", completed_work, total_work, "", genome_scope or "all", stored_rows)
    return DirectMatrixBuildSummary(
        output_file=output_file,
        sample_count=len(sample_ids),
        stored_rows=stored_rows,
        storage_layout=CURRENT_MATRIX_HDF5_SPARSE_LAYOUT if sparse else CURRENT_MATRIX_HDF5_LAYOUT,
        genome_scope=genome_scope or "all",
    )


def append_matrix_hdf5_from_duckdb(
    conn: duckdb.DuckDBPyConnection,
    *,
    sample_ids: list[str],
    matrix_hdf5_file: Path,
    memory_limit_gb: float = 16.0,
    duckdb_threads: int = 1,
    progress_callback: BuildProgressCallback | None = None,
) -> DirectMatrixAppendSummary:
    """Append samples directly from MetaTrawl DuckDB rows into an existing HDF5 matrix."""
    if not sample_ids:
        raise ValueError("No samples were provided for matrix append.")
    _configure_duckdb_for_matrix(conn, memory_limit_gb=memory_limit_gb, threads=duckdb_threads)
    matrix_hdf5_file = Path(matrix_hdf5_file).expanduser().resolve()
    if not matrix_hdf5_file.exists():
        raise FileNotFoundError(f"Matrix store file does not exist: {matrix_hdf5_file}")

    h5py = _import_h5py()
    memory_limit_bytes = _memory_limit_bytes(memory_limit_gb)
    stored_rows = 0
    completed_work = 0
    with h5py.File(str(matrix_hdf5_file), "r+") as handle:
        metadata = {str(k): str(v) for k, v in handle["metadata"].attrs.items()}
        count_dtype = str(metadata.get("count_dtype", "uint16"))
        layout = str(metadata.get("layout", CURRENT_MATRIX_HDF5_LAYOUT))
        genome_scope = str(metadata.get("genome_scope", "all"))
        sparse = layout == CURRENT_MATRIX_HDF5_SPARSE_LAYOUT
        genomes = _read_genomes(handle["genomes"])
        genome_scaffolds = _read_genome_scaffolds(handle["genome_scaffolds"])
        existing_samples = _read_samples(handle["samples"])
        existing_names = {name for _idx, name in existing_samples}
        duplicates = sorted(set(sample_ids) & existing_names)
        if duplicates:
            raise ValueError("Cannot append samples already present in matrix store: " + ", ".join(duplicates))
        observed = _observed_genomes_for_samples(
            conn,
            sample_ids=sample_ids,
            genome=None if genome_scope == "all" else genome_scope,
        )
        allowed = {spec.genome for spec in genomes}
        missing = observed - allowed
        if missing:
            raise ValueError("Profiles contain genomes missing from existing matrix store: " + ", ".join(sorted(missing)))

        total_sample_count = len(existing_samples) + len(sample_ids)
        existing_sample_count = len(existing_samples)
        sample_rows = [(existing_sample_count + offset, sample_id) for offset, sample_id in enumerate(sample_ids)]
        _resize_sample_datasets(handle["samples"], existing_samples=existing_samples, sample_rows=sample_rows)
        total_work = len(sample_rows) * len(genomes)
        _emit_progress(progress_callback, "start", completed_work, total_work, "", genome_scope, stored_rows)
        scaffolds_by_genome_idx = _group_scaffolds_by_genome(genome_scaffolds)
        sparse_states: dict[int, tuple[object, object, int]] = {}
        for spec in genomes:
            _check_matrix_memory(spec, count_dtype=count_dtype, memory_limit_bytes=memory_limit_bytes)
            node = handle[_hdf5_matrix_dataset_path(spec.genome_idx)]
            if sparse:
                indptr = node["indptr"]
                indices = node["indices"]
                indptr.resize((total_sample_count + 1,))
                sparse_states[spec.genome_idx] = (indptr, indices, int(indptr[existing_sample_count]))
            else:
                node.resize((total_sample_count, spec.matrix_length, 4))

        for sample_row, sample_id in sample_rows:
            for spec in genomes:
                _emit_progress(progress_callback, "processing", completed_work, total_work, sample_id, spec.genome, stored_rows)
                matrix = _load_duckdb_sample_genome_matrix(
                    conn,
                    sample_id=sample_id,
                    genome_spec=spec,
                    genome_offsets=scaffolds_by_genome_idx[spec.genome_idx],
                    count_dtype=count_dtype,
                )
                if sparse:
                    indptr_ds, indices_ds, current_nnz = sparse_states[spec.genome_idx]
                    current_nnz = _append_sparse_hdf5_matrix_row(
                        indptr_dataset=indptr_ds,
                        indices_dataset=indices_ds,
                        sample_row=sample_row,
                        flat_indices=_dense_matrix_to_sparse_flat_indices(matrix),
                        current_nnz=current_nnz,
                    )
                    sparse_states[spec.genome_idx] = (indptr_ds, indices_ds, current_nnz)
                else:
                    handle[_hdf5_matrix_dataset_path(spec.genome_idx)][sample_row, :, :] = matrix
                stored_rows += 1
                completed_work += 1
                _emit_progress(progress_callback, "advance", completed_work, total_work, sample_id, spec.genome, stored_rows)
        _emit_progress(progress_callback, "done", completed_work, total_work, "", genome_scope, stored_rows)
    return DirectMatrixAppendSummary(
        output_file=matrix_hdf5_file,
        appended_sample_count=len(sample_ids),
        total_sample_count=total_sample_count,
        stored_rows=stored_rows,
    )


def _load_duckdb_sample_genome_matrix(
    conn: duckdb.DuckDBPyConnection,
    *,
    sample_id: str,
    genome_spec: GenomeSpec,
    genome_offsets: list[GenomeScaffoldOffset],
    count_dtype: str,
) -> np.ndarray:
    np_dtype = COUNT_DTYPES[count_dtype]
    matrix = np.zeros((genome_spec.matrix_length, 4), dtype=np_dtype)
    rows = conn.execute(
        """
        SELECT chrom, pos, A, T, C, G
        FROM profile_positions
        WHERE sample_id = ? AND genome = ?
          AND CAST(A AS BIGINT) + CAST(T AS BIGINT) + CAST(C AS BIGINT) + CAST(G AS BIGINT) >= ?
        """,
        [sample_id, genome_spec.genome, MATRIX_BUILD_MIN_COV],
    ).fetchall()
    if not rows:
        return matrix
    offsets_by_chrom = {offset.chrom: offset for offset in genome_offsets}
    for chrom, pos, a, t, c, g in rows:
        offset = offsets_by_chrom.get(str(chrom))
        if offset is None:
            raise ValueError(f"Profile sample {sample_id} has scaffold {chrom} outside matrix contract for genome {genome_spec.genome}.")
        pos_int = int(pos)
        if pos_int < offset.min_pos or pos_int > offset.max_pos:
            raise ValueError(
                f"Profile sample {sample_id} has out-of-range position {genome_spec.genome}:{chrom}:{pos_int}; "
                f"matrix range is {offset.min_pos}-{offset.max_pos}."
            )
        axis_pos = pos_int - offset.index_base + offset.axis_start
        matrix[axis_pos, 0] = 1 if int(a) > 0 else 0
        matrix[axis_pos, 1] = 1 if int(t) > 0 else 0
        matrix[axis_pos, 2] = 1 if int(c) > 0 else 0
        matrix[axis_pos, 3] = 1 if int(g) > 0 else 0
    return matrix


def _initialize_hdf5_store(
    handle,
    *,
    h5py_module,
    metadata_rows: dict[str, str],
    sample_rows: list[tuple[int, str]],
    genomes: list[GenomeSpec],
    genome_scaffolds: list[GenomeScaffoldOffset],
    contract_genomes: list[GenomeSpec],
    contract_genome_scaffolds: list[GenomeScaffoldOffset],
    gene_ranges: list[GeneRangeSpec],
    count_dtype: str,
    export_batch_mb: float,
    sparse: bool,
) -> None:
    handle.attrs["zipstrain_hdf5_version"] = MATRIX_HDF5_FILE_VERSION
    metadata = handle.create_group("metadata")
    for key, value in metadata_rows.items():
        metadata.attrs[str(key)] = str(value)
    sample_count = len(sample_rows)
    sample_chunk_len = _matrix_hdf5_sample_axis_chunk_length(sample_count)
    samples_group = handle.create_group("samples")
    samples_group.create_dataset(
        "sample_idx",
        data=np.asarray([sample_idx for sample_idx, _sample_name in sample_rows], dtype=np.int64),
        chunks=(sample_chunk_len,),
        maxshape=(None,),
    )
    _write_hdf5_string_dataset(
        samples_group,
        "sample_name",
        [sample_name for _sample_idx, sample_name in sample_rows],
        h5py_module=h5py_module,
        chunks=(sample_chunk_len,),
        maxshape=(None,),
    )
    _write_hdf5_genomes_group(handle, "genomes", genomes, h5py_module=h5py_module)
    _write_hdf5_genome_scaffolds_group(handle, "genome_scaffolds", genome_scaffolds, h5py_module=h5py_module)
    _write_hdf5_genomes_group(handle, MATRIX_HDF5_CONTRACT_GENOMES_GROUP, contract_genomes, h5py_module=h5py_module)
    _write_hdf5_genome_scaffolds_group(
        handle,
        MATRIX_HDF5_CONTRACT_SCAFFOLDS_GROUP,
        contract_genome_scaffolds,
        h5py_module=h5py_module,
    )
    if gene_ranges:
        _write_hdf5_gene_ranges_group(handle, "genes", gene_ranges, h5py_module=h5py_module)
    matrices_group = handle.create_group("matrices")
    for spec in genomes:
        if sparse:
            _create_sparse_hdf5_genome_store(
                matrices_group,
                genome_idx=spec.genome_idx,
                sample_count=sample_count,
                matrix_length=spec.matrix_length,
            )
        else:
            chunk_samples = _matrix_hdf5_chunk_sample_count(
                sample_count=sample_count,
                matrix_length=spec.matrix_length,
                dtype_name=count_dtype,
                target_batch_mb=export_batch_mb,
            )
            matrices_group.create_dataset(
                str(spec.genome_idx),
                shape=(sample_count, spec.matrix_length, 4),
                dtype=COUNT_DTYPES[count_dtype],
                chunks=(chunk_samples, spec.matrix_length, 4),
                maxshape=(None, spec.matrix_length, 4),
                fillvalue=0,
            )


def _metadata_rows(
    *,
    genome: str,
    count_dtype: str,
    memory_limit_gb: float,
    export_batch_mb: float,
    sparse: bool,
    bed_file: Path,
    stb_file: Path,
    gene_range_table: Path | None,
    has_gene_ranges: bool,
) -> dict[str, str]:
    rows = {
        "profiles_dir": "metatrawl_duckdb",
        "profile_format": "metatrawl_duckdb_profile_positions",
        "genome_scope": genome,
        "count_dtype": count_dtype,
        "layout": CURRENT_MATRIX_HDF5_SPARSE_LAYOUT if sparse else CURRENT_MATRIX_HDF5_LAYOUT,
        "matrix_value_semantics": FILTERED_PRESENCE_MATRIX_VALUE_SEMANTICS,
        "coverage_filter_min_cov": str(MATRIX_BUILD_MIN_COV),
        "memory_limit_gb": str(memory_limit_gb),
        "export_batch_mb": str(export_batch_mb),
        "separator_rows_between_scaffolds": "1",
        "input_format": "hdf5",
        "has_gene_ranges": "1" if has_gene_ranges else "0",
        "bed_file": str(Path(bed_file).resolve()),
        "stb_file": str(Path(stb_file).resolve()),
    }
    if gene_range_table is not None:
        rows["gene_range_table"] = str(Path(gene_range_table).resolve())
    return rows


def _collect_scaffold_specs_from_bed_and_stb(*, bed_file: Path, stb_file: Path, genome: str | None) -> list[ScaffoldSpec]:
    bed_path = Path(bed_file)
    if not bed_path.is_file():
        raise FileNotFoundError(f"BED file does not exist: {bed_path}")
    stb_path = Path(stb_file)
    if not stb_path.is_file():
        raise FileNotFoundError(f"STB file does not exist: {stb_path}")
    bed_spans = (
        pl.scan_csv(bed_path, separator="\t", has_header=False)
        .select(
            pl.col("column_1").cast(pl.Utf8).alias("chrom"),
            pl.col("column_2").cast(pl.Int64).alias("start"),
            pl.col("column_3").cast(pl.Int64).alias("end"),
        )
        .group_by("chrom")
        .agg(pl.col("start").min().alias("min_start"), pl.col("end").max().alias("max_end"))
        .collect(engine="streaming")
    )
    stb_mapping = _read_stb_mapping(stb_path)
    if genome is not None:
        stb_mapping = stb_mapping.filter(pl.col("genome") == genome)
    if stb_mapping.height == 0:
        detail = genome if genome is not None else "all"
        raise ValueError(f"No STB rows found for genome scope: {detail}")
    contract_df = stb_mapping.join(bed_spans, left_on="scaffold", right_on="chrom", how="left")
    missing = contract_df.filter(pl.col("max_end").is_null() | pl.col("min_start").is_null()).get_column("scaffold").to_list()
    if missing:
        raise ValueError("BED file is missing scaffold intervals for: " + ", ".join(sorted(str(x) for x in missing)))
    specs: list[ScaffoldSpec] = []
    for scaffold_idx, row in enumerate(contract_df.sort(["genome", "scaffold"]).iter_rows(named=True)):
        min_start = int(row["min_start"])
        max_end = int(row["max_end"])
        index_base = min_start + 1
        vector_length = max_end - min_start
        if vector_length <= 0:
            raise ValueError(f"Invalid BED span for scaffold {row['scaffold']}: start={min_start}, end={max_end}")
        specs.append(
            ScaffoldSpec(
                scaffold_idx=scaffold_idx,
                genome=str(row["genome"]),
                chrom=str(row["scaffold"]),
                index_base=index_base,
                vector_length=vector_length,
                min_pos=index_base,
                max_pos=max_end,
            )
        )
    return specs


def _read_stb_mapping(stb_file: Path) -> pl.DataFrame:
    mapping = (
        pl.scan_csv(stb_file, separator="\t", has_header=False)
        .select(pl.col("column_1").cast(pl.Utf8).alias("scaffold"), pl.col("column_2").cast(pl.Utf8).alias("genome"))
        .collect(engine="streaming")
    )
    duplicate_scaffolds = mapping.group_by("scaffold").len().filter(pl.col("len") > 1).get_column("scaffold").to_list()
    if duplicate_scaffolds:
        raise ValueError("STB file must map each scaffold to exactly one genome. Duplicate scaffold mappings found for: " + ", ".join(sorted(str(x) for x in duplicate_scaffolds)))
    return mapping.unique(["scaffold", "genome"])


def _build_genome_specs(scaffolds: list[ScaffoldSpec]) -> tuple[list[GenomeSpec], list[GenomeScaffoldOffset]]:
    grouped: dict[str, list[ScaffoldSpec]] = {}
    for spec in scaffolds:
        grouped.setdefault(spec.genome, []).append(spec)
    genomes: list[GenomeSpec] = []
    offsets: list[GenomeScaffoldOffset] = []
    for genome_idx, genome_name in enumerate(sorted(grouped)):
        axis_cursor = 0
        true_length = 0
        genome_scaffolds = sorted(grouped[genome_name], key=lambda item: (item.chrom, item.index_base))
        for ordinal, spec in enumerate(genome_scaffolds):
            axis_start = axis_cursor
            axis_end = axis_start + spec.vector_length - 1
            offsets.append(
                GenomeScaffoldOffset(
                    genome_idx=genome_idx,
                    scaffold_ordinal=ordinal,
                    genome=genome_name,
                    chrom=spec.chrom,
                    axis_start=axis_start,
                    axis_end=axis_end,
                    index_base=spec.index_base,
                    vector_length=spec.vector_length,
                    min_pos=spec.min_pos,
                    max_pos=spec.max_pos,
                )
            )
            axis_cursor = axis_end + 1
            true_length += spec.vector_length
            if ordinal < len(genome_scaffolds) - 1:
                axis_cursor += 1
        genomes.append(GenomeSpec(genome_idx=genome_idx, genome=genome_name, matrix_length=axis_cursor, true_length=true_length, scaffold_count=len(genome_scaffolds)))
    return genomes, offsets


def _subset_contract_to_genomes(*, contract_genomes: list[GenomeSpec], contract_genome_scaffolds: list[GenomeScaffoldOffset], genome_names: set[str]) -> tuple[list[GenomeSpec], list[GenomeScaffoldOffset]]:
    return (
        [spec for spec in contract_genomes if spec.genome in genome_names],
        [offset for offset in contract_genome_scaffolds if offset.genome in genome_names],
    )


def _collect_gene_range_specs(*, gene_range_table: Path, genome_scaffolds: list[GenomeScaffoldOffset]) -> list[GeneRangeSpec]:
    if not gene_range_table.is_file():
        raise FileNotFoundError(f"Gene range table does not exist: {gene_range_table}")
    offsets_by_chrom = {offset.chrom: offset for offset in genome_scaffolds}
    frame = (
        pl.scan_csv(gene_range_table, has_header=False, separator="\t")
        .rename({"column_1": "gene", "column_2": "scaffold", "column_3": "start", "column_4": "end"})
        .select(pl.col("gene").cast(pl.Utf8), pl.col("scaffold").cast(pl.Utf8), pl.col("start").cast(pl.Int64), pl.col("end").cast(pl.Int64))
        .collect(engine="streaming")
    )
    specs: list[GeneRangeSpec] = []
    for row in frame.iter_rows(named=True):
        chrom = str(row["scaffold"])
        offset = offsets_by_chrom.get(chrom)
        if offset is None:
            continue
        start = int(row["start"])
        end = int(row["end"])
        if end < start:
            raise ValueError(f"Gene {row['gene']} on scaffold {chrom} has end < start ({end} < {start}).")
        if start < offset.min_pos or end > offset.max_pos:
            raise ValueError(f"Gene {row['gene']} on scaffold {chrom} falls outside the matrix coordinate range: gene={start}-{end}, matrix={offset.min_pos}-{offset.max_pos}")
        specs.append(GeneRangeSpec(-1, str(row["gene"]), offset.genome_idx, offset.genome, chrom, offset.axis_start + start - offset.index_base, offset.axis_start + end - offset.index_base))
    sorted_specs = sorted(specs, key=lambda item: (item.genome_idx, item.axis_start, item.axis_end, item.gene))
    return [GeneRangeSpec(idx, spec.gene, spec.genome_idx, spec.genome, spec.chrom, spec.axis_start, spec.axis_end) for idx, spec in enumerate(sorted_specs)]


def _observed_genomes_for_samples(conn: duckdb.DuckDBPyConnection, *, sample_ids: list[str], genome: str | None) -> set[str]:
    if not sample_ids:
        return set()
    conditions = ["sample_id IN (SELECT unnest(?))"]
    params: list[object] = [sample_ids]
    if genome is not None:
        conditions.append("genome = ?")
        params.append(genome)
    rows = conn.execute(f"SELECT DISTINCT genome FROM profile_positions WHERE {' AND '.join(conditions)}", params).fetchall()
    return {str(row[0]) for row in rows}


def _configure_duckdb_for_matrix(conn: duckdb.DuckDBPyConnection, *, memory_limit_gb: float | None, threads: int) -> None:
    conn.execute(f"SET threads = {max(1, int(threads))}")
    if memory_limit_gb is not None:
        conn.execute(f"SET memory_limit = {_sql_literal(f'{memory_limit_gb}GB')}")


def _write_hdf5_string_dataset(group, name: str, values: list[str], *, h5py_module, maxshape=None, chunks=None) -> None:
    kwargs = {"data": np.asarray(values, dtype=object), "dtype": h5py_module.string_dtype(encoding="utf-8")}
    if maxshape is not None:
        kwargs["maxshape"] = maxshape
    if chunks is not None:
        kwargs["chunks"] = chunks
    group.create_dataset(name, **kwargs)


def _read_hdf5_string_dataset(dataset) -> list[str]:
    if hasattr(dataset, "asstr"):
        return [str(value) for value in dataset.asstr()[...].tolist()]
    raw = dataset[...]
    return [value.decode("utf-8") if isinstance(value, (bytes, bytearray, np.bytes_)) else str(value) for value in raw.tolist()]


def _write_hdf5_genomes_group(handle, group_name: str, genomes: list[GenomeSpec], *, h5py_module) -> None:
    group = handle.create_group(group_name)
    group.create_dataset("genome_idx", data=np.asarray([spec.genome_idx for spec in genomes], dtype=np.int64))
    _write_hdf5_string_dataset(group, "genome", [spec.genome for spec in genomes], h5py_module=h5py_module)
    group.create_dataset("matrix_length", data=np.asarray([spec.matrix_length for spec in genomes], dtype=np.int64))
    group.create_dataset("true_length", data=np.asarray([spec.true_length for spec in genomes], dtype=np.int64))
    group.create_dataset("scaffold_count", data=np.asarray([spec.scaffold_count for spec in genomes], dtype=np.int64))


def _write_hdf5_genome_scaffolds_group(handle, group_name: str, genome_scaffolds: list[GenomeScaffoldOffset], *, h5py_module) -> None:
    group = handle.create_group(group_name)
    group.create_dataset("genome_idx", data=np.asarray([spec.genome_idx for spec in genome_scaffolds], dtype=np.int64))
    group.create_dataset("scaffold_ordinal", data=np.asarray([spec.scaffold_ordinal for spec in genome_scaffolds], dtype=np.int64))
    _write_hdf5_string_dataset(group, "genome", [spec.genome for spec in genome_scaffolds], h5py_module=h5py_module)
    _write_hdf5_string_dataset(group, "chrom", [spec.chrom for spec in genome_scaffolds], h5py_module=h5py_module)
    for field in ("axis_start", "axis_end", "index_base", "vector_length", "min_pos", "max_pos"):
        group.create_dataset(field, data=np.asarray([getattr(spec, field) for spec in genome_scaffolds], dtype=np.int64))


def _write_hdf5_gene_ranges_group(handle, group_name: str, gene_ranges: list[GeneRangeSpec], *, h5py_module) -> None:
    group = handle.create_group(group_name)
    group.create_dataset("gene_idx", data=np.asarray([spec.gene_idx for spec in gene_ranges], dtype=np.int64))
    group.create_dataset("genome_idx", data=np.asarray([spec.genome_idx for spec in gene_ranges], dtype=np.int64))
    _write_hdf5_string_dataset(group, "genome", [spec.genome for spec in gene_ranges], h5py_module=h5py_module)
    _write_hdf5_string_dataset(group, "chrom", [spec.chrom for spec in gene_ranges], h5py_module=h5py_module)
    _write_hdf5_string_dataset(group, "gene", [spec.gene for spec in gene_ranges], h5py_module=h5py_module)
    group.create_dataset("axis_start", data=np.asarray([spec.axis_start for spec in gene_ranges], dtype=np.int64))
    group.create_dataset("axis_end", data=np.asarray([spec.axis_end for spec in gene_ranges], dtype=np.int64))


def _read_genomes(group) -> list[GenomeSpec]:
    return [
        GenomeSpec(int(idx), str(genome), int(matrix_length), int(true_length), int(scaffold_count))
        for idx, genome, matrix_length, true_length, scaffold_count in zip(
            np.asarray(group["genome_idx"][...], dtype=np.int64).tolist(),
            _read_hdf5_string_dataset(group["genome"]),
            np.asarray(group["matrix_length"][...], dtype=np.int64).tolist(),
            np.asarray(group["true_length"][...], dtype=np.int64).tolist(),
            np.asarray(group["scaffold_count"][...], dtype=np.int64).tolist(),
        )
    ]


def _read_genome_scaffolds(group) -> list[GenomeScaffoldOffset]:
    return [
        GenomeScaffoldOffset(int(genome_idx), int(scaffold_ordinal), str(genome), str(chrom), int(axis_start), int(axis_end), int(index_base), int(vector_length), int(min_pos), int(max_pos))
        for genome_idx, scaffold_ordinal, genome, chrom, axis_start, axis_end, index_base, vector_length, min_pos, max_pos in zip(
            np.asarray(group["genome_idx"][...], dtype=np.int64).tolist(),
            np.asarray(group["scaffold_ordinal"][...], dtype=np.int64).tolist(),
            _read_hdf5_string_dataset(group["genome"]),
            _read_hdf5_string_dataset(group["chrom"]),
            np.asarray(group["axis_start"][...], dtype=np.int64).tolist(),
            np.asarray(group["axis_end"][...], dtype=np.int64).tolist(),
            np.asarray(group["index_base"][...], dtype=np.int64).tolist(),
            np.asarray(group["vector_length"][...], dtype=np.int64).tolist(),
            np.asarray(group["min_pos"][...], dtype=np.int64).tolist(),
            np.asarray(group["max_pos"][...], dtype=np.int64).tolist(),
        )
    ]


def _read_samples(group) -> list[tuple[int, str]]:
    return list(zip(np.asarray(group["sample_idx"][...], dtype=np.int64).astype(int).tolist(), _read_hdf5_string_dataset(group["sample_name"])))


def _resize_sample_datasets(group, *, existing_samples: list[tuple[int, str]], sample_rows: list[tuple[int, str]]) -> None:
    total = len(existing_samples) + len(sample_rows)
    group["sample_idx"].resize((total,))
    group["sample_name"].resize((total,))
    group["sample_idx"][len(existing_samples):total] = np.asarray([idx for idx, _name in sample_rows], dtype=np.int64)
    group["sample_name"][len(existing_samples):total] = np.asarray([name for _idx, name in sample_rows], dtype=object)


def _hdf5_matrix_dataset_path(genome_idx: int) -> str:
    return f"matrices/{int(genome_idx)}"


def _matrix_hdf5_sample_axis_chunk_length(sample_count: int) -> int:
    return max(1, min(int(sample_count), 1024))


def _matrix_hdf5_chunk_sample_count(*, sample_count: int, matrix_length: int, dtype_name: str, target_batch_mb: float) -> int:
    dtype = COUNT_DTYPES[dtype_name]
    per_sample_bytes = max(1, matrix_length * 4 * np.dtype(dtype).itemsize)
    target_batch_bytes = int(target_batch_mb * (1024 ** 2))
    return max(1, min(sample_count, int(target_batch_bytes // per_sample_bytes) or 1))


def _matrix_hdf5_sparse_indices_chunk_length(matrix_length: int) -> int:
    return max(1024, min(max(matrix_length * 4, 1), 1_048_576))


def _create_sparse_hdf5_genome_store(matrices_group, *, genome_idx: int, sample_count: int, matrix_length: int):
    group = matrices_group.create_group(str(genome_idx))
    group.create_dataset("indptr", shape=(sample_count + 1,), dtype=np.int64, chunks=(_matrix_hdf5_sample_axis_chunk_length(sample_count + 1),), maxshape=(None,), fillvalue=0)
    group.create_dataset("indices", shape=(0,), dtype=np.int64, chunks=(_matrix_hdf5_sparse_indices_chunk_length(matrix_length),), maxshape=(None,), fillvalue=0)
    return group


def _initial_sparse_states(handle, genomes: list[GenomeSpec]) -> dict[int, tuple[object, object, int]]:
    states = {}
    for spec in genomes:
        group = handle[_hdf5_matrix_dataset_path(spec.genome_idx)]
        states[spec.genome_idx] = (group["indptr"], group["indices"], 0)
    return states


def _dense_matrix_to_sparse_flat_indices(matrix: np.ndarray) -> np.ndarray:
    return np.flatnonzero(matrix.reshape(-1)).astype(np.int64, copy=False)


def _append_sparse_hdf5_matrix_row(*, indptr_dataset, indices_dataset, sample_row: int, flat_indices: np.ndarray, current_nnz: int) -> int:
    flat_indices = np.asarray(flat_indices, dtype=np.int64)
    next_nnz = current_nnz + int(flat_indices.size)
    if flat_indices.size > 0:
        indices_dataset.resize((next_nnz,))
        indices_dataset[current_nnz:next_nnz] = flat_indices
    indptr_dataset[int(sample_row) + 1] = next_nnz
    return next_nnz


def _group_scaffolds_by_genome(genome_scaffolds: list[GenomeScaffoldOffset]) -> dict[int, list[GenomeScaffoldOffset]]:
    grouped: dict[int, list[GenomeScaffoldOffset]] = {}
    for offset in genome_scaffolds:
        grouped.setdefault(offset.genome_idx, []).append(offset)
    return grouped


def _memory_limit_bytes(memory_limit_gb: float) -> int:
    return int(float(memory_limit_gb) * (1024 ** 3))


def _check_matrix_memory(spec: GenomeSpec, *, count_dtype: str, memory_limit_bytes: int) -> None:
    estimated = spec.matrix_length * 4 * np.dtype(COUNT_DTYPES[count_dtype]).itemsize
    if estimated > memory_limit_bytes:
        raise MemoryError(f"Genome {spec.genome} matrix requires about {estimated} bytes, exceeding configured memory limit {memory_limit_bytes} bytes.")


def _emit_progress(callback: BuildProgressCallback | None, phase: str, completed: int, total: int, sample_name: str, genome: str, stored_rows: int) -> None:
    if callback is None:
        return
    callback({"phase": phase, "completed": completed, "total": total, "sample_name": sample_name, "genome": genome, "scaffold": "", "stored_rows": stored_rows})


def _import_h5py():
    try:
        import h5py  # type: ignore
    except ImportError as exc:
        raise RuntimeError("HDF5 matrix operations require h5py. Install MetaTrawl with matrix dependencies.") from exc
    return h5py


def _sql_literal(value: object) -> str:
    return "'" + str(value).replace("'", "''") + "'"

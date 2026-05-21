"""Workflow adapters around ZipStrain matrix operations."""

from __future__ import annotations

from pathlib import Path
import shutil
from tempfile import TemporaryDirectory

from metatrawl import db


def build_matrix_from_registry(
    conn,
    *,
    output_file: Path,
    genome: str,
    overwrite: bool = False,
    matrix_id: str | None = None,
    count_dtype: str = "uint16",
    memory_limit_gb: float = 16.0,
    export_batch_mb: float = 128.0,
) -> db.MatrixStore:
    """Build a ZipStrain matrix store from all complete profiles in the registry."""
    profile_paths = db.complete_profile_paths(conn)
    if not profile_paths:
        raise ValueError("No complete profiles are available for matrix build.")

    output_file = Path(output_file)
    if output_file.exists() and not overwrite:
        raise FileExistsError(f"Matrix output already exists: {output_file}")
    output_file.parent.mkdir(parents=True, exist_ok=True)

    chosen_matrix_id = matrix_id or output_file.stem
    if db.get_matrix_store(conn, chosen_matrix_id) is not None and not overwrite:
        raise ValueError(f"Matrix ID already exists: {chosen_matrix_id}")

    from zipstrain import matrix_pairs as mp

    with TemporaryDirectory(prefix="metatrawl_profiles_") as tmp_dir:
        profile_dir = Path(tmp_dir)
        for idx, profile_path in enumerate(profile_paths, start=1):
            link_name = _unique_profile_name(profile_path, idx)
            _link_or_copy(profile_path, profile_dir / link_name)

        mp.build_matrix_hdf5(
            profile_dir=profile_dir,
            output_file=output_file,
            genome=genome,
            count_dtype=count_dtype,
            memory_limit_gb=memory_limit_gb,
            export_batch_mb=export_batch_mb,
        )

    return db.register_matrix_store(
        conn,
        matrix_id=chosen_matrix_id,
        genome=genome,
        matrix_file=output_file,
        profile_count=len(profile_paths),
        overwrite=overwrite,
    )


def compare_matrix_from_registry(
    conn,
    *,
    matrix_id: str,
    output_file: Path,
    calculate: str = "all",
    genome: str = "all",
    backend: str = "numpy",
    memory_limit_gb: float = 16.0,
) -> str:
    """Run ZipStrain matrix compare for a registered matrix store."""
    store = db.get_matrix_store(conn, matrix_id)
    if store is None:
        raise ValueError(f"Unknown matrix ID: {matrix_id}")

    output_file = Path(output_file)
    output_file.parent.mkdir(parents=True, exist_ok=True)

    from zipstrain import matrix_pairs as mp

    mp.matrix_compare(
        matrix_db_file=store.matrix_file,
        output_file=output_file,
        genome=genome,
        memory_limit_gb=memory_limit_gb,
        backend=backend,
        calculate=calculate,
    )

    compare_id = output_file.stem
    db.register_matrix_compare(
        conn,
        compare_id=compare_id,
        matrix_id=matrix_id,
        compare_db_file=output_file,
        calculate=calculate,
    )
    return compare_id


def _unique_profile_name(profile_path: Path, idx: int) -> str:
    suffix = profile_path.suffix or ".parquet"
    stem = profile_path.stem or f"profile_{idx}"
    return f"{idx:06d}_{stem}{suffix}"


def _link_or_copy(src: Path, dst: Path) -> None:
    try:
        dst.symlink_to(src.resolve())
    except OSError:
        shutil.copy2(src, dst)

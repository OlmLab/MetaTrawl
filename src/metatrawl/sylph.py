"""Read Sylph results without assuming diagnostic columns are numeric."""

from pathlib import Path

import polars as pl

from metatrawl.logging import WorkflowLogger


def read_sylph_table(path: Path) -> pl.DataFrame:
    """Read one sample, preserving mixed types and removing exact repeated headers.

    Sylph legitimately emits HIGH/LOW in Eff_lambda. Numeric validation belongs
    to consumers of abundance/coverage, not to unrelated diagnostic columns.
    """
    path = Path(path)
    if path.suffix.lower() == ".parquet":
        table = pl.read_parquet(path)
    else:
        separator = "\t" if path.suffix.lower() in {".tsv", ".txt", ".tab"} else ","
        table = pl.read_csv(path, separator=separator, infer_schema=False)
    if table.is_empty():
        return table
    repeated_header = pl.all_horizontal([
        pl.col(name).cast(pl.String).eq(name).fill_null(False) for name in table.columns
    ])
    cleaned = table.filter(~repeated_header)
    skipped = table.height - cleaned.height
    if skipped:
        WorkflowLogger().emit(step="sylph-table", status="repeated-header-skipped", path=path, rows=skipped)
    if "Sample_file" in cleaned.columns and cleaned["Sample_file"].drop_nulls().n_unique() > 1:
        raise ValueError(
            f"Sylph table contains multiple Sample_file values: {path}. "
            "MetaTrawl expects one Sylph profile per sample. Regenerate the abundance file "
            "with a single FASTQ or with R1/R2 passed together as -1 and -2."
        )
    return cleaned


def numeric_column(table: pl.DataFrame, column: str) -> pl.Series:
    """Reject malformed required numbers rather than silently dropping records."""
    try:
        values = table[column].cast(pl.Float64, strict=True)
    except pl.exceptions.PolarsError as exc:
        raise ValueError(f"Invalid numeric values in Sylph column {column!r}: {exc}") from exc
    if values.is_nan().any() or values.is_infinite().any() or (values < 0).any():
        raise ValueError(f"Sylph column {column!r} must contain finite, non-negative numbers.")
    return values

"""Compact, fixed-threshold allele-mask profile storage.

ZipStrain still produces its normal profile parquet. MetaTrawl converts that
temporary full profile into covered-position and allele-set masks during import.
The stored masks use ZipStrain's HDF5 bit order: A=1, T=2, C=4, G=8.
"""

from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
from pathlib import Path
from typing import Iterable, Iterator

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
import zstandard as zstd


PROFILE_STORAGE_FULL = "full"
PROFILE_STORAGE_ALLELE_MASK = "allele-mask"
PROFILE_STORAGE_MODES = (PROFILE_STORAGE_FULL, PROFILE_STORAGE_ALLELE_MASK)
ALLELE_MASK_FORMAT_VERSION = 1
ALLELE_MASK_BIT_ORDER = "ATCG"
ALLELE_MASK_CODEC = "zstd"
PROFILE_BATCH_ROWS = 262_144

_COMPRESSOR = zstd.ZstdCompressor(level=6, write_checksum=True)
_DECOMPRESSOR = zstd.ZstdDecompressor()
_REFERENCE_BITS = np.zeros(256, dtype=np.uint8)
for _base, _value in ((b"A", 1), (b"T", 2), (b"C", 4), (b"G", 8)):
    _REFERENCE_BITS[_base[0]] = _value
    _REFERENCE_BITS[_base.lower()[0]] = _value


@dataclass(frozen=True)
class AlleleMaskImportSummary:
    """Counts written while compacting one full profile."""

    blocks: int
    covered_positions: int
    source_rows: int


@dataclass(frozen=True)
class ReferenceSegment:
    """One shared reference scaffold."""

    segment_id: int
    genome: str
    chrom: str
    start_pos: int
    span: int
    masks: np.ndarray


@dataclass(frozen=True)
class DecodedProfileBlock:
    """One decoded sample/scaffold allele-mask block."""

    chrom: str
    start_pos: int
    masks: np.ndarray


def store_profile_parquet(
    conn,
    *,
    sample_id: str,
    profile_file: Path,
    min_cov: int,
    cache_dir: Path | None,
) -> AlleleMaskImportSummary:
    """Encode a normal ZipStrain profile parquet into allele-mask tables."""
    profile_file = Path(profile_file)
    if profile_file.suffix.lower() != ".parquet":
        raise ValueError("Allele-mask storage requires a ZipStrain profile parquet.")
    parquet_file = pq.ParquetFile(profile_file)
    required = {"chrom", "genome", "pos", "A", "C", "G", "T"}
    missing = required - set(parquet_file.schema_arrow.names)
    if missing:
        raise ValueError(
            "profile_file missing required columns for allele-mask storage: "
            + ", ".join(sorted(missing))
        )
    batches = parquet_file.iter_batches(
        batch_size=PROFILE_BATCH_ROWS,
        columns=["chrom", "genome", "pos", "A", "C", "G", "T"],
    )
    return store_profile_batches(
        conn,
        sample_id=sample_id,
        batches=batches,
        min_cov=min_cov,
        cache_dir=cache_dir,
    )


def store_profile_batches(
    conn,
    *,
    sample_id: str,
    batches: Iterable[pa.RecordBatch],
    min_cov: int,
    cache_dir: Path | None,
) -> AlleleMaskImportSummary:
    """Stream ordered profile batches into compact sample/scaffold blocks."""
    if min_cov < 1:
        raise ValueError("Allele-mask min_cov must be >= 1.")

    current_key: tuple[str, str] | None = None
    current_segment: ReferenceSegment | None = None
    current_presence: np.ndarray | None = None
    current_deviation: np.ndarray | None = None
    current_last_pos: int | None = None
    completed_keys: set[tuple[str, str]] = set()
    ensured_genomes: set[str] = set()
    segment_cache: dict[tuple[str, str], ReferenceSegment] = {}
    block_count = 0
    covered_total = 0
    source_total = 0

    def finish_current() -> None:
        nonlocal block_count, covered_total
        if current_segment is None or current_presence is None or current_deviation is None:
            return
        covered = int(current_presence.sum())
        if covered == 0:
            return
        presence_packed = pack_presence(current_presence)
        deviation_packed = pack_nibbles(current_deviation)
        payload_hash = sha256(presence_packed + deviation_packed).hexdigest()
        conn.execute(
            """
            INSERT INTO allele_mask_profile_blocks
              (sample_id, segment_id, presence, deviation, covered_positions, payload_hash)
            VALUES (?, ?, ?, ?, ?, ?)
            """,
            [
                sample_id,
                current_segment.segment_id,
                _compress(presence_packed),
                _compress(deviation_packed),
                covered,
                payload_hash,
            ],
        )
        block_count += 1
        covered_total += covered

    for batch in batches:
        if batch.num_rows == 0:
            continue
        source_total += batch.num_rows
        names = set(batch.schema.names)
        missing = {"chrom", "genome", "pos", "A", "C", "G", "T"} - names
        if missing:
            raise ValueError(
                "Profile batch missing required columns: " + ", ".join(sorted(missing))
            )
        chroms = np.asarray(batch.column("chrom").to_pylist(), dtype=object)
        genomes = np.asarray(batch.column("genome").to_pylist(), dtype=object)
        positions = _numeric_column(batch, "pos")
        counts = {
            base: _numeric_column(batch, base)
            for base in ("A", "C", "G", "T")
        }
        if np.any(positions < 1):
            raise ValueError("Profile positions must use positive, one-based coordinates.")
        if any(np.any(values < 0) for values in counts.values()):
            raise ValueError("Profile allele counts cannot be negative.")

        boundaries = np.flatnonzero(
            (chroms[1:] != chroms[:-1]) | (genomes[1:] != genomes[:-1])
        ) + 1
        starts = np.concatenate(([0], boundaries))
        stops = np.concatenate((boundaries, [batch.num_rows]))
        for start, stop in zip(starts.tolist(), stops.tolist()):
            genome = str(genomes[start])
            chrom = str(chroms[start])
            key = (genome, chrom)
            if key != current_key:
                if current_key is not None:
                    finish_current()
                    completed_keys.add(current_key)
                if key in completed_keys:
                    raise ValueError(
                        "ZipStrain profile rows must be grouped by scaffold for "
                        f"allele-mask import; scaffold repeated: {genome}/{chrom}"
                    )
                current_key = key
                current_last_pos = None
                if genome in {"", "NA", "None", "null"}:
                    current_segment = None
                    current_presence = None
                    current_deviation = None
                else:
                    if genome not in ensured_genomes:
                        ensure_reference_genome(
                            conn,
                            genome=genome,
                            cache_dir=cache_dir,
                        )
                        ensured_genomes.add(genome)
                    current_segment = segment_cache.get(key)
                    if current_segment is None:
                        current_segment = load_reference_segment(
                            conn,
                            genome=genome,
                            chrom=chrom,
                        )
                        segment_cache[key] = current_segment
                    current_presence = np.zeros(current_segment.span, dtype=np.bool_)
                    current_deviation = np.zeros(current_segment.span, dtype=np.uint8)

            if current_segment is None:
                continue
            group_positions = positions[start:stop]
            if (
                (current_last_pos is not None and int(group_positions[0]) <= current_last_pos)
                or np.any(np.diff(group_positions) <= 0)
            ):
                raise ValueError(
                    "ZipStrain profile positions must be strictly increasing within "
                    f"each scaffold; invalid order in {genome}/{chrom}."
                )
            current_last_pos = int(group_positions[-1])
            offsets = group_positions - current_segment.start_pos
            if np.any(offsets < 0) or np.any(offsets >= current_segment.span):
                bad_pos = int(group_positions[(offsets < 0) | (offsets >= current_segment.span)][0])
                raise ValueError(
                    f"Profile position {genome}/{chrom}:{bad_pos} is outside cached "
                    f"reference range {current_segment.start_pos}-"
                    f"{current_segment.start_pos + current_segment.span - 1}."
                )
            coverage = (
                counts["A"][start:stop]
                + counts["T"][start:stop]
                + counts["C"][start:stop]
                + counts["G"][start:stop]
            )
            covered = coverage >= min_cov
            if not np.any(covered):
                continue
            allele_masks = (
                (counts["A"][start:stop] > 0).astype(np.uint8)
                | ((counts["T"][start:stop] > 0).astype(np.uint8) << 1)
                | ((counts["C"][start:stop] > 0).astype(np.uint8) << 2)
                | ((counts["G"][start:stop] > 0).astype(np.uint8) << 3)
            )
            covered_offsets = offsets[covered].astype(np.int64, copy=False)
            covered_masks = allele_masks[covered]
            if np.any(covered_masks == 0):
                raise ValueError(
                    f"Covered profile positions produced an empty allele mask in {genome}/{chrom}."
                )
            current_presence[covered_offsets] = True
            current_deviation[covered_offsets] = (
                covered_masks ^ current_segment.masks[covered_offsets]
            )

    finish_current()
    return AlleleMaskImportSummary(
        blocks=block_count,
        covered_positions=covered_total,
        source_rows=source_total,
    )


def ensure_reference_genome(conn, *, genome: str, cache_dir: Path | None) -> None:
    """Store all reference scaffolds for a genome once, validating existing rows."""
    existing = int(
        conn.execute(
            "SELECT count(*) FROM allele_mask_reference_segments WHERE genome = ?",
            [genome],
        ).fetchone()[0]
    )
    if existing:
        return
    fasta = _resolve_genome_fasta(conn, genome=genome, cache_dir=cache_dir)
    records = list(_read_fasta(fasta))
    if not records:
        raise ValueError(f"Cached genome FASTA contains no sequences: {fasta}")
    next_segment_id = int(
        conn.execute(
            "SELECT COALESCE(max(segment_id), 0) + 1 FROM allele_mask_reference_segments"
        ).fetchone()[0]
    )
    for ordinal, (chrom, sequence) in enumerate(records):
        masks = encode_reference_sequence(sequence)
        packed = pack_nibbles(masks)
        conn.execute(
            """
            INSERT INTO allele_mask_reference_segments
              (segment_id, genome, chrom, segment_ordinal, start_pos, span,
               reference_mask, reference_hash)
            VALUES (?, ?, ?, ?, 1, ?, ?, ?)
            """,
            [
                next_segment_id + ordinal,
                genome,
                chrom,
                ordinal,
                int(masks.size),
                _compress(packed),
                sha256(packed).hexdigest(),
            ],
        )


def load_reference_segment(conn, *, genome: str, chrom: str) -> ReferenceSegment:
    """Load and verify one shared reference scaffold."""
    row = conn.execute(
        """
        SELECT segment_id, start_pos, span, reference_mask, reference_hash
        FROM allele_mask_reference_segments
        WHERE genome = ? AND chrom = ?
        """,
        [genome, chrom],
    ).fetchone()
    if row is None:
        raise ValueError(
            f"Cached reference for genome {genome} does not contain scaffold {chrom}."
        )
    span = int(row[2])
    packed = _decompress(row[3], expected_size=(span + 1) // 2)
    if sha256(packed).hexdigest() != str(row[4]):
        raise ValueError(f"Reference checksum mismatch for {genome}/{chrom}.")
    return ReferenceSegment(
        segment_id=int(row[0]),
        genome=genome,
        chrom=chrom,
        start_pos=int(row[1]),
        span=span,
        masks=unpack_nibbles(packed, span),
    )


def iter_decoded_profile_blocks(
    conn,
    *,
    sample_id: str,
    genome: str,
) -> Iterator[DecodedProfileBlock]:
    """Yield verified allele-mask blocks for one sample and genome."""
    rows = conn.execute(
        """
        SELECT r.chrom, r.start_pos, r.span, r.reference_mask, r.reference_hash,
               p.presence, p.deviation, p.covered_positions, p.payload_hash
        FROM allele_mask_profile_blocks p
        JOIN allele_mask_reference_segments r USING (segment_id)
        WHERE p.sample_id = ? AND r.genome = ?
        ORDER BY r.segment_ordinal, r.chrom
        """,
        [sample_id, genome],
    ).fetchall()
    for (
        chrom,
        start_pos,
        span_value,
        reference_blob,
        reference_hash,
        presence_blob,
        deviation_blob,
        covered_positions,
        payload_hash,
    ) in rows:
        span = int(span_value)
        reference_packed = _decompress(
            reference_blob,
            expected_size=(span + 1) // 2,
        )
        if sha256(reference_packed).hexdigest() != str(reference_hash):
            raise ValueError(f"Reference checksum mismatch for {genome}/{chrom}.")
        presence_packed = _decompress(
            presence_blob,
            expected_size=(span + 7) // 8,
        )
        deviation_packed = _decompress(
            deviation_blob,
            expected_size=(span + 1) // 2,
        )
        if sha256(presence_packed + deviation_packed).hexdigest() != str(payload_hash):
            raise ValueError(
                f"Allele-mask payload checksum mismatch for {sample_id}/{genome}/{chrom}."
            )
        presence = unpack_presence(presence_packed, span)
        if int(presence.sum()) != int(covered_positions):
            raise ValueError(
                f"Allele-mask covered-position count mismatch for "
                f"{sample_id}/{genome}/{chrom}."
            )
        reference = unpack_nibbles(reference_packed, span)
        deviation = unpack_nibbles(deviation_packed, span)
        masks = np.zeros(span, dtype=np.uint8)
        masks[presence] = reference[presence] ^ deviation[presence]
        if np.any(masks[presence] == 0):
            raise ValueError(
                f"Allele-mask payload decoded an empty covered allele set for "
                f"{sample_id}/{genome}/{chrom}."
            )
        yield DecodedProfileBlock(
            chrom=str(chrom),
            start_pos=int(start_pos),
            masks=masks,
        )


def encode_reference_sequence(sequence: str) -> np.ndarray:
    """Encode FASTA bases using the matrix A/T/C/G bit order."""
    raw = np.frombuffer(sequence.encode("ascii"), dtype=np.uint8)
    return _REFERENCE_BITS[raw]


def pack_presence(values: np.ndarray) -> bytes:
    """Pack one Boolean value per position using little-endian bit order."""
    return np.packbits(
        np.asarray(values, dtype=np.bool_),
        bitorder="little",
    ).tobytes()


def unpack_presence(payload: bytes, length: int) -> np.ndarray:
    """Unpack a presence bitmap to exactly ``length`` positions."""
    values = np.unpackbits(
        np.frombuffer(payload, dtype=np.uint8),
        bitorder="little",
    )
    return values[:length].astype(np.bool_, copy=False)


def pack_nibbles(values: np.ndarray) -> bytes:
    """Pack two 4-bit masks per byte."""
    masks = np.asarray(values, dtype=np.uint8)
    if np.any(masks > 15):
        raise ValueError("Allele masks must fit in four bits.")
    packed = np.zeros((masks.size + 1) // 2, dtype=np.uint8)
    packed[:] = masks[0::2]
    if masks.size > 1:
        packed[: masks.size // 2] |= masks[1::2] << 4
    return packed.tobytes()


def unpack_nibbles(payload: bytes, length: int) -> np.ndarray:
    """Unpack exactly ``length`` four-bit masks."""
    packed = np.frombuffer(payload, dtype=np.uint8)
    values = np.empty(packed.size * 2, dtype=np.uint8)
    values[0::2] = packed & 0x0F
    values[1::2] = packed >> 4
    return values[:length]


def _numeric_column(batch: pa.RecordBatch, name: str) -> np.ndarray:
    column = batch.column(name)
    if column.null_count:
        raise ValueError(f"Profile column {name} contains null values.")
    values = column.to_numpy(zero_copy_only=False)
    if np.issubdtype(values.dtype, np.floating):
        if np.any(~np.isfinite(values)) or np.any(values != np.floor(values)):
            raise ValueError(f"Profile column {name} must contain integer values.")
    return values.astype(np.int64, copy=False)


def _compress(payload: bytes) -> bytes:
    return _COMPRESSOR.compress(payload)


def _decompress(payload: bytes, *, expected_size: int) -> bytes:
    try:
        value = _DECOMPRESSOR.decompress(payload, max_output_size=expected_size)
    except zstd.ZstdError as exc:
        raise ValueError("Invalid zstd payload in allele-mask storage.") from exc
    if len(value) != expected_size:
        raise ValueError(
            f"Invalid allele-mask payload size: expected {expected_size}, got {len(value)}."
        )
    return value


def _resolve_genome_fasta(conn, *, genome: str, cache_dir: Path | None) -> Path:
    row = conn.execute(
        "SELECT genome_fasta FROM cache_genomes WHERE accession = ?",
        [genome],
    ).fetchone()
    candidates: list[Path] = []
    if row is not None and row[0]:
        candidates.append(Path(str(row[0])))
    if cache_dir is not None:
        candidates.append(
            Path(cache_dir) / "genomes" / f"{_safe_name(genome)}.fna"
        )
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    detail = ", ".join(str(path) for path in candidates) or "no cache directory provided"
    raise FileNotFoundError(
        f"Allele-mask import needs the cached FASTA for genome {genome}; checked: {detail}"
    )


def _read_fasta(path: Path) -> Iterator[tuple[str, str]]:
    chrom: str | None = None
    chunks: list[str] = []
    with Path(path).open() as handle:
        for line in handle:
            if line.startswith(">"):
                if chrom is not None:
                    yield chrom, "".join(chunks)
                chrom = line[1:].strip().split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
    if chrom is not None:
        yield chrom, "".join(chunks)


def _safe_name(value: str) -> str:
    return "".join(char if char.isalnum() or char in "._-" else "_" for char in value)

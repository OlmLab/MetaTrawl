"""Small, transactional database and profiling provenance records.

This module never reads profile_positions. Unknown historical settings remain
unknown; acknowledging a mismatch never changes an existing contract.
"""
from __future__ import annotations

from dataclasses import asdict
from functools import lru_cache
import hashlib
from importlib.metadata import PackageNotFoundError, version
import json
from pathlib import Path
import time
import uuid

SCHEMA_VERSION = 1
MANIFEST_VERSION = 1


def canonical(value) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


def fingerprint(value) -> str:
    return hashlib.sha256(canonical(value).encode()).hexdigest()


def package_versions() -> dict:
    result = {}
    for name in ("metatrawl", "zipstrain", "duckdb"):
        try:
            result[name] = version(name)
        except PackageNotFoundError:
            result[name] = "unknown"
    return result


def table_exists(conn, name: str) -> bool:
    return conn.execute("SELECT 1 FROM information_schema.tables WHERE table_schema='main' AND table_name=?", [name]).fetchone() is not None


def check_schema(conn) -> None:
    """Reject unsupported schemas before any DDL or writable migration."""
    if not table_exists(conn, "database_metadata"):
        return
    rows = conn.execute("SELECT schema_version FROM database_metadata").fetchall()
    if len(rows) != 1 or rows[0][0] != SCHEMA_VERSION:
        raise ValueError(f"Unsupported MetaTrawl database schema: {rows!r}; supported={SCHEMA_VERSION}. Use a compatible MetaTrawl release.")


def migrate(conn) -> None:
    """Install provenance metadata once, touching only the sample registry."""
    if table_exists(conn, "database_metadata"):
        check_schema(conn)
        return
    conn.execute("BEGIN TRANSACTION")
    try:
        conn.execute('''
            CREATE TABLE database_metadata (
                id INTEGER PRIMARY KEY CHECK (id=1), database_uuid VARCHAR NOT NULL,
                schema_version INTEGER NOT NULL, created_at DOUBLE NOT NULL,
                created_by VARCHAR NOT NULL, origin VARCHAR NOT NULL,
                accepted_contract_id VARCHAR);
            CREATE TABLE profiling_contracts (
                contract_id VARCHAR PRIMARY KEY, parameters_json VARCHAR NOT NULL,
                created_at DOUBLE NOT NULL);
            CREATE TABLE workflow_runs (
                run_id VARCHAR PRIMARY KEY, started_at DOUBLE NOT NULL, finished_at DOUBLE,
                status VARCHAR NOT NULL, contract_id VARCHAR, details_json VARCHAR NOT NULL);
            CREATE TABLE sample_provenance (
                sample_id VARCHAR PRIMARY KEY, contract_id VARCHAR, provenance_status VARCHAR NOT NULL,
                workflow_run_id VARCHAR, manifest_json VARCHAR NOT NULL, updated_at DOUBLE NOT NULL);
            CREATE TABLE profile_references (
                genome VARCHAR NOT NULL, reference_id VARCHAR NOT NULL,
                details_json VARCHAR NOT NULL, PRIMARY KEY(genome, reference_id));
            CREATE TABLE profiling_contract_assets (
                contract_id VARCHAR NOT NULL, asset_kind VARCHAR NOT NULL,
                asset_hash VARCHAR NOT NULL, PRIMARY KEY(contract_id, asset_kind));
            CREATE TABLE provenance_events (
                event_id VARCHAR PRIMARY KEY, created_at DOUBLE NOT NULL,
                sample_id VARCHAR, workflow_run_id VARCHAR, kind VARCHAR NOT NULL, details_json VARCHAR NOT NULL);
        ''')
        legacy = conn.execute("SELECT EXISTS(SELECT 1 FROM samples)").fetchone()[0]
        now = time.time()
        conn.execute("INSERT INTO database_metadata VALUES (1, ?, ?, ?, ?, ?, NULL)",
                     [str(uuid.uuid4()), SCHEMA_VERSION, now, package_versions()["metatrawl"], "legacy-adopted" if legacy else "new"])
        conn.execute("INSERT INTO sample_provenance SELECT sample_id, NULL, 'legacy-unknown', NULL, '{}', ? FROM samples", [now])
        conn.execute("INSERT INTO metatrawl_migrations VALUES ('provenance_schema_v1', ?)", [now])
        conn.execute("COMMIT")
    except Exception:
        conn.execute("ROLLBACK")
        raise


@lru_cache(maxsize=1)
def zipstrain_defaults() -> dict:
    # Resolve defaults from the installed producer, not today's MetaTrawl defaults.
    from zipstrain import profile, utils
    return {
        "min_read_ani": profile.PROFILE_MIN_READ_ANI_DEFAULT,
        "null_model": {"error_rate": utils.NULL_MODEL_ERROR_RATE_DEFAULT,
                       "max_total_reads": utils.NULL_MODEL_MAX_COVERAGE_DEFAULT,
                       "p_threshold": utils.NULL_MODEL_P_THRESHOLD_DEFAULT,
                       "model_type": "poisson"},
    }


def profiling_contract(config) -> dict:
    values = asdict(config)
    defaults = zipstrain_defaults()
    if values["min_read_ani"] is None:
        values["min_read_ani"] = defaults["min_read_ani"]
    return {"contract_version": 1, "profile": values,
            "null_model": defaults["null_model"],
            "alignment": {"aligner": "bowtie2", "preset": "default-sensitive-end-to-end",
                          "samtools_exclude_flags": 4},
            "single_end_read_inclusion": "all-mapped"}


def validate_contract(contract: dict) -> None:
    if not isinstance(contract, dict) or contract.get("contract_version") != 1:
        raise ValueError("Unsupported or missing profiling contract_version (expected 1).")
    required = {"min_mapq", "min_baseq", "min_freq", "min_read_ani", "read_inclusion"}
    if not isinstance(contract.get("profile"), dict) or not required <= contract["profile"].keys():
        raise ValueError("Profiling contract must contain all effective profile settings.")
    if not all(key in contract for key in ("null_model", "alignment", "single_end_read_inclusion")):
        raise ValueError("Profiling contract is missing null_model, alignment, or single_end_read_inclusion.")
    profile = contract["profile"]
    for key in ("min_mapq", "min_baseq"):
        if type(profile[key]) is not int or profile[key] < 0:
            raise ValueError(f"Invalid profiling contract {key}.")
    for key in ("min_freq", "min_read_ani"):
        if type(profile[key]) not in (int, float) or not 0 <= profile[key] <= 1:
            raise ValueError(f"Invalid profiling contract {key}.")
    if profile["read_inclusion"] not in ("paired", "proper-pairs", "all-mapped"):
        raise ValueError("Invalid profiling contract read_inclusion.")
    model = contract["null_model"]
    if not isinstance(model, dict) or not {"error_rate", "max_total_reads", "p_threshold", "model_type"} <= model.keys():
        raise ValueError("Incomplete null-model settings in profiling contract.")
    if not isinstance(contract["alignment"], dict):
        raise ValueError("Invalid alignment settings in profiling contract.")
    canonical(contract)


def differences(expected, actual, prefix="") -> list[str]:
    if isinstance(expected, dict) and isinstance(actual, dict):
        result = []
        for key in sorted(expected.keys() | actual.keys()):
            result.extend(differences(expected.get(key), actual.get(key), f"{prefix}.{key}" if prefix else key))
        return result
    return [] if expected == actual else [f"{prefix}: database={expected!r} incoming={actual!r}"]


def compatibility_issues(conn, contract: dict | None) -> list[str]:
    issues = []
    if contract is None:
        issues.append("incoming profile has unknown historical settings (no provenance manifest)")
    else:
        validate_contract(contract)
    baseline = conn.execute("SELECT accepted_contract_id FROM database_metadata WHERE id=1").fetchone()[0]
    if baseline and contract is not None:
        expected = json.loads(conn.execute("SELECT parameters_json FROM profiling_contracts WHERE contract_id=?", [baseline]).fetchone()[0])
        issues.extend(differences(expected, contract))
    if contract is not None and conn.execute("SELECT 1 FROM sample_provenance WHERE contract_id IS NULL LIMIT 1").fetchone():
        issues.append("existing samples have legacy-unknown profiling settings")
    return issues


def require_compatible(issues: list[str], allow: bool) -> None:
    if issues and not allow:
        raise ValueError("Profiling compatibility cannot be established:\n  " + "\n  ".join(issues) +
                         "\nUse --allow-incompatible-profiles to acknowledge these differences; they will be recorded.")


def event(conn, *, kind, details, sample_id=None, run_id=None) -> None:
    conn.execute("INSERT INTO provenance_events VALUES (?, ?, ?, ?, ?, ?)",
                 [str(uuid.uuid4()), time.time(), sample_id, run_id, kind, canonical(details)])


def register_contract(conn, contract) -> str:
    validate_contract(contract)
    identity = fingerprint(contract)
    conn.execute("INSERT OR IGNORE INTO profiling_contracts VALUES (?, ?, ?)", [identity, canonical(contract), time.time()])
    conn.execute("UPDATE database_metadata SET accepted_contract_id=? WHERE id=1 AND accepted_contract_id IS NULL", [identity])
    return identity


def start_run(conn, *, contract, details, allow=False) -> str:
    """Preflight and record one workflow invocation before expensive work starts."""
    conn.execute("BEGIN TRANSACTION")
    try:
        issues = compatibility_issues(conn, contract)
        require_compatible(issues, allow)
        identity = register_contract(conn, contract)
        run_id = str(uuid.uuid4())
        conn.execute("INSERT INTO workflow_runs VALUES (?, ?, NULL, 'running', ?, ?)",
                     [run_id, time.time(), identity, canonical(details)])
        if issues:
            event(conn, kind="compatibility-override", details=issues, run_id=run_id)
        conn.execute("COMMIT")
        return run_id
    except Exception:
        conn.execute("ROLLBACK")
        raise


def finish_run(conn, run_id, status) -> None:
    conn.execute("UPDATE workflow_runs SET status=?, finished_at=? WHERE run_id=?", [status, time.time(), run_id])


def manifest_path(profile_file: Path) -> Path:
    return Path(str(profile_file) + ".provenance.json")


def read_manifest(profile_file: Path, sample_id: str, explicit_path: Path | None = None) -> dict | None:
    path = explicit_path or manifest_path(profile_file)
    if not path.exists():
        if explicit_path:
            raise FileNotFoundError(f"Provenance manifest does not exist: {path}")
        return None
    value = json.loads(path.read_text())
    if not isinstance(value, dict) or value.get("manifest_version") != MANIFEST_VERSION or value.get("sample_id") != sample_id:
        raise ValueError(f"Invalid provenance manifest version or sample identity: {path}")
    validate_contract(value.get("contract"))
    if value.get("contract_id") != fingerprint(value["contract"]):
        raise ValueError(f"Provenance contract hash mismatch: {path}")
    if value.get("status") not in ("recorded", "user-declared"):
        raise ValueError(f"Invalid provenance status: {path}")
    if not isinstance(value.get("references"), dict) or not isinstance(value.get("details"), dict):
        raise ValueError(f"Invalid provenance references/details: {path}")
    for reference in value["references"].values():
        if not isinstance(reference, dict) or not {"sequence_hash", "annotation_hash", "scaffolds"} <= reference.keys():
            raise ValueError(f"Incomplete genome reference identity: {path}")
    return value


def record_import(conn, *, sample_id, manifest, allow=False, run_id=None, storage_mode="full") -> None:
    """Called inside the same transaction as profile rows and completion markers."""
    contract = manifest["contract"] if manifest else None
    issues = compatibility_issues(conn, contract)
    references = manifest.get("references", {}) if manifest else {}
    null_hash = manifest.get("details", {}).get("null_model_sha256") if manifest else None
    identity = fingerprint(contract) if contract else None
    if null_hash:
        expected = conn.execute("SELECT asset_hash FROM profiling_contract_assets WHERE contract_id=? AND asset_kind='null-model'", [identity]).fetchone()
        if expected and expected[0] != null_hash:
            issues.append(f"null-model content differs for identical profiling parameters: database={expected[0]} incoming={null_hash}")
    for genome, reference in references.items():
        previous = conn.execute("SELECT details_json FROM profile_references WHERE genome=?", [genome]).fetchall()
        changed = [json.loads(row[0]) for row in previous if json.loads(row[0]) != reference]
        if changed:
            if storage_mode == "allele-mask" and any(old.get("sequence_hash") != reference.get("sequence_hash") for old in changed):
                raise ValueError(f"Cannot change allele-mask reference sequence for {genome}, even with an override.")
            issues.append(f"reference sequence/scaffolds or gene annotation differ for genome {genome}")
    require_compatible(issues, allow)
    identity = register_contract(conn, contract) if contract else None
    if null_hash:
        conn.execute("INSERT OR IGNORE INTO profiling_contract_assets VALUES (?, 'null-model', ?)", [identity, null_hash])
    old = conn.execute("SELECT manifest_json FROM sample_provenance WHERE sample_id=?", [sample_id]).fetchone()
    if old:
        event(conn, kind="sample-provenance-replaced", details=json.loads(old[0]), sample_id=sample_id, run_id=run_id)
    conn.execute("INSERT OR REPLACE INTO sample_provenance VALUES (?, ?, ?, ?, ?, ?)",
                 [sample_id, identity, manifest["status"] if manifest else "legacy-unknown", run_id,
                  canonical(manifest or {}), time.time()])
    for genome, reference in references.items():
        conn.execute("INSERT OR IGNORE INTO profile_references VALUES (?, ?, ?)", [genome, fingerprint(reference), canonical(reference)])
    if issues:
        event(conn, kind="compatibility-override", details=issues, sample_id=sample_id, run_id=run_id)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def fasta_records(path, *, full_header=False):
    name = None
    digest = None
    length = 0
    with Path(path).open() as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, length, digest.hexdigest()
                name = line[1:] if full_header else line[1:].split()[0]
                digest, length = hashlib.sha256(), 0
            else:
                if name is None:
                    raise ValueError(f"Invalid FASTA: {path}")
                sequence = line.upper().encode("ascii")
                digest.update(sequence)
                length += len(sequence)
    if name is not None:
        yield name, length, digest.hexdigest()


def reference_catalog(reference) -> dict:
    """Fingerprint per genome, independent of sample reference concatenation order."""
    mapping = {}
    for line in Path(reference.stb_file).read_text().splitlines():
        if line.strip():
            chrom, genome = line.split()[:2]
            if chrom in mapping:
                raise ValueError(f"Duplicate scaffold in STB: {chrom}")
            mapping[chrom] = genome
    genomes = {}
    for chrom, length, digest in fasta_records(reference.reference_fasta):
        if chrom not in mapping:
            raise ValueError(f"Reference scaffold missing from STB: {chrom}")
        entry = genomes.setdefault(mapping[chrom], {"scaffolds": [], "genes": []})
        entry["scaffolds"].append([chrom, length, digest])
    for header, length, digest in fasta_records(reference.gene_fasta, full_header=True):
        gene = header.split()[0]
        chrom = gene.rsplit("_", 1)[0]
        if chrom not in mapping:
            raise ValueError(f"Gene scaffold missing from STB: {gene}")
        genomes[mapping[chrom]]["genes"].append([gene, length, digest, header])
    return {genome: {"sequence_hash": fingerprint(sorted(entry["scaffolds"])),
                     "annotation_hash": fingerprint(sorted(entry["genes"])),
                     "scaffolds": sorted(entry["scaffolds"])} for genome, entry in genomes.items()}


def write_manifest(*, profile_file, sample_id, contract, references, details, status="recorded") -> None:
    payload = {"manifest_version": MANIFEST_VERSION, "sample_id": sample_id,
               "status": status, "contract": contract, "contract_id": fingerprint(contract),
               "references": references, "details": details, "created_at": time.time()}
    target = manifest_path(profile_file)
    temporary = target.with_name(f".{target.name}.{uuid.uuid4().hex}.tmp")
    temporary.write_text(canonical(payload) + "\n")
    temporary.replace(target)


def null_model_options() -> list[str]:
    values = zipstrain_defaults()["null_model"]
    return ["--error-rate", str(values["error_rate"]), "--max-total-reads", str(values["max_total_reads"]),
            "--p-threshold", str(values["p_threshold"]), "--model-type", values["model_type"]]


def describe(conn) -> dict:
    if not table_exists(conn, "database_metadata"):
        return {"schema_version": "legacy-unversioned", "provenance": "legacy-unknown"}
    check_schema(conn)
    names = [item[0] for item in conn.execute("SELECT * FROM database_metadata LIMIT 0").description]
    result = dict(zip(names, conn.execute("SELECT * FROM database_metadata").fetchone()))
    if table_exists(conn, "profile_storage"):
        storage = conn.execute("SELECT mode, format_version, min_cov, codec FROM profile_storage WHERE id=1").fetchone()
        if storage:
            result["profile_storage"] = dict(zip(("mode", "format_version", "min_cov", "codec"), storage))
    result["sample_provenance"] = dict(conn.execute("SELECT provenance_status, count(*) FROM sample_provenance GROUP BY provenance_status").fetchall())
    result["contracts"] = [{"contract_id": row[0], "parameters": json.loads(row[1])} for row in conn.execute("SELECT contract_id, parameters_json FROM profiling_contracts ORDER BY contract_id").fetchall()]
    return result


def operational_config(config) -> dict:
    """Keep execution settings without persisting potentially secret env values."""
    result = asdict(config)
    for stage in result["stages"].values():
        stage["environment"] = {key: "<set>" for key in stage.get("environment", {})}
    return result


def declare_legacy(conn, *, contract, samples, reason) -> int:
    """Attach user-declared settings without claiming historical verification."""
    if not reason.strip():
        raise ValueError("A reason is required for a historical declaration.")
    validate_contract(contract)
    conn.execute("BEGIN TRANSACTION")
    try:
        rows = conn.execute("SELECT sample_id FROM sample_provenance WHERE contract_id IS NULL ORDER BY sample_id").fetchall()
        eligible = {row[0] for row in rows}
        selected = eligible if samples is None else set(samples)
        if not selected <= eligible:
            raise ValueError("Declarations can only target existing legacy-unknown samples.")
        identity = register_contract(conn, contract)
        for sample in sorted(selected):
            manifest = {"status": "user-declared", "contract": contract, "reason": reason}
            conn.execute("UPDATE sample_provenance SET contract_id=?, provenance_status='user-declared', manifest_json=?, updated_at=? WHERE sample_id=?",
                         [identity, canonical(manifest), time.time(), sample])
            event(conn, kind="legacy-contract-declared", details=manifest, sample_id=sample)
        conn.execute("COMMIT")
        return len(selected)
    except Exception:
        conn.execute("ROLLBACK")
        raise

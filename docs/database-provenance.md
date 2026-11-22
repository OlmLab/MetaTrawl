# Database versions and profiling provenance

MetaTrawl records the database schema separately from its Python package version.
A database receives a permanent UUID when it is created. Renaming, copying, or
moving the file preserves the UUID; it is an identity, not a checksum. Two
copies can subsequently diverge.

## Older databases are read-only

**MetaTrawl 2 writes only to databases that MetaTrawl 2 created.** Every older
database stays readable, but any command that would write to it stops with an
error before changing the file:

| Database | Created by | MetaTrawl 2 | Keep writing with |
| --- | --- | --- | --- |
| Schema v2, origin `new` | MetaTrawl 2 | read and write | MetaTrawl 2 |
| Schema v1 | MetaTrawl 1.x | read-only | MetaTrawl 1.x |
| Unversioned (no `database_metadata` table) | MetaTrawl 0.x | read-only | MetaTrawl 0.x or 1.x |
| Schema v2, origin `legacy-adopted` | MetaTrawl 2.0.x opening a 0.x file | read-only | none; see below |

Read-only access covers `database info`, `database history`, `runs list`,
`profiles list`, `profiles remaining`, `status`, and the [Python API](python-api.md),
including pair analyses. Older files have no recorded profiling settings, so
pass `allow_unknown_provenance=True` to `pair(...)` to acknowledge that:

```python
from metatrawl import open_database
pair = open_database("old_project.duckdb").pair(
    "SRR1", "SRR2", genome="GCF_000005845", allow_unknown_provenance=True)
```

To keep profiling an older project with its original settings, use the release
or container that created it. To profile with ZipStrain 2, start a new database.
There is no in-place upgrade: an older database holds ZipStrain 1 profiles made
with different defaults (for example `min_freq`), and relabelling the file would
not change how those profiles were made.

MetaTrawl 2.0.0 and 2.0.1 adopted unversioned files as schema v2 on their first
write. Such files are marked `legacy-adopted`, MetaTrawl 1.x refuses to write
to them, and MetaTrawl 2.1 opens them read-only. A backup taken before the
adoption still works with MetaTrawl 0.x and 1.x.

A file is created as a new v2 database only when it is missing or empty. A
newer, unsupported schema is rejected before modification, even for reads.
Overrides cannot bypass any of these checks.

## Inspect the database

```sh
metatrawl database info --db metatrawl.duckdb
metatrawl database history --db metatrawl.duckdb --limit 50
```

`info` is read-only and prints the UUID, schema version, storage representation,
profiling contracts, and provenance-status counts. `history` lists explicit
overrides and sample-provenance replacements.

```python
from metatrawl import open_database
metadata = open_database("metatrawl.duckdb").metadata()
```

The `workflow_runs`, `profiling_contracts`, `sample_provenance`,
`profile_references`, and `provenance_events` tables are ordinary queryable
DuckDB tables. Workflow history records resolved operational settings, start/end
times, and completion status. Environment values are redacted. An abruptly killed
process can leave a run marked `running`; sample imports remain transactional.

## New profiling runs

`sync-profile` resolves effective settings before starting workers. The first
accepted contract establishes the project's baseline. Later runs with different
quality thresholds, allele-frequency filters, read policies, alignment settings,
or null-model settings are rejected with a field-by-field explanation.

Worker counts, RAM, queues, retries, and Slurm settings are operational settings,
not scientific incompatibilities. The intentional single-ended `all-mapped`
override is recorded per sample without conflicting with the requested paired
policy. Different sample-specific Sylph reference sets are allowed. Reference
sequences, scaffold names/lengths, and Prodigal gene identities/coordinates are
fingerprinted per genome, so the same genome name cannot silently change meaning.
Null-model file identities are also checked when available.

Optional `[profile.null_model]` fields and `--null-model-*` CLI overrides are
resolved against the installed ZipStrain defaults before the contract is
recorded. Unset overrides preserve existing default contract identities.
Model changes are scientific changes, unlike RAM or Slurm walltime changes.

After profiling, outputs include a sidecar such as:

```text
SRR123.profile.parquet.provenance.json
```

Keep it with the bundle when transferring outputs. It contains the settings used
to generate that sample, its reference identities, read layout, command and source
information. Installed Python package versions are recorded from the coordinator;
this does not independently attest the software installed inside a different
Slurm stage environment. Package-version differences alone do not reject data;
contract and reference differences do. Future known semantic incompatibilities
require updating the contract policy.

V2 contracts explicitly include the profile generation, `collect_codons`, and
`allele_integration`. The importer also authenticates ZipStrain's v2 catalog,
reference-codon hash, and optional read-phased counts. Generation changes are
scientific changes, not an operational backend switch.

The importer records the manifest in the same transaction as profile/stat rows
and completion status. Cleanup deletes the sidecar only after a successful commit.
Missing historical sidecars remain unknown rather than receiving current settings.
`profiles import` automatically finds the sidecar, or accepts `--provenance-file`.
`profiles add` accepts an optional `provenance_file` column in its bundle CSV.

## Explicit overrides

To add data when historical compatibility is unknown, or deliberately accept a
scientific mismatch, add this flag to `sync-profile`, `profiles import`, or
`profiles add`:

```text
--allow-incompatible-profiles
```

The differences are recorded permanently. The baseline contract and historical
sample settings are not replaced. This is acknowledgment, not a claim that mixed
profiles produce scientifically equivalent results. Malformed manifests,
unsupported formats, and changed reference sequences in allele-mask storage
cannot be forced through.

Old incomplete scratch checkpoints also require acknowledgment when their
settings are unknown. Such continued samples are labeled `user-declared`, not
`recorded`. A checkpoint with a known conflicting contract must be regenerated
in a new scratch directory, even with an import override. Published old bundles
are imported as-is with unknown provenance; resuming does not relabel them.

## Contract template

```sh
metatrawl database contract-template --workflow-config project.toml > contract.json
```

prints the profiling contract that a workflow configuration resolves to on the
**current installation**. Use it to compare settings between configurations;
it is not a reconstruction of the defaults an older release used.

## Scope

This first version enforces compatibility at profiling preflight and DuckDB
import. Existing matrix append, comparison resume, and genome-view checks remain
unchanged. Propagating
per-sample provenance into HDF5/compare artifacts and enforcing those identities
is a separate downstream step. A UUID alone does not establish compatibility.

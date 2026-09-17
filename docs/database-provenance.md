# Database versions and profiling provenance

MetaTrawl records the database schema separately from its Python package version.
A database receives a permanent UUID when it is created or adopted. Renaming,
copying, or moving the file preserves the UUID; it is an identity, not a checksum.
Two copies can subsequently diverge.

## Existing databases

Existing databases remain readable without migration. The first writable open
adds small provenance tables and labels registered samples `legacy-unknown`.
It does not infer old parameters from current defaults, scan profile positions,
or convert the existing profile count columns. Existing profile, stats, matrix,
and comparison data do not need to be rebuilt for this metadata adoption.

You can adopt explicitly:

```sh
metatrawl database migrate --db metatrawl.duckdb
```

Use a backup before upgrading a valuable database. Older MetaTrawl releases do
not know about the new guards; do not use them to write to an adopted database.
A newer, unsupported schema is rejected before modification. Overrides cannot
bypass this protection. The metadata creation timestamp on an adopted database
is the adoption time, not an inferred historical creation date.

## Inspect the database

```sh
metatrawl database info --db metatrawl.duckdb
metatrawl database history --db metatrawl.duckdb --limit 50
```

`info` is read-only and prints the UUID, schema version, storage representation,
profiling contracts, and provenance-status counts. `history` lists explicit
overrides, sample-provenance replacements, and historical declarations.

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

## Declare known historical settings

If you have the original configuration, you can attach it to legacy samples:

```sh
metatrawl database contract-template --workflow-config historical.toml > contract.json
# Inspect every field against the historical software and configuration first.
metatrawl database declare-contract \
  --db metatrawl.duckdb \
  --contract-file contract.json \
  --reason "Checked against archived workflow configuration"
```

The template resolves defaults from the **current installation**; it is not a
reconstruction of old defaults. Correct the JSON from historical evidence before
using it. Repeat `--sample-id` to restrict the declaration; otherwise it applies
to all `legacy-unknown` samples. Declared settings stay labeled `user-declared`.
Reference sequences are not reconstructed by this declaration.

## Scope

This first version enforces compatibility at profiling preflight and DuckDB
import. Existing matrix append, comparison resume, and genome-view checks remain
unchanged; database adoption alone does not invalidate old artifacts. Propagating
per-sample provenance into HDF5/compare artifacts and enforcing those identities
is a separate downstream step. A UUID alone does not establish compatibility.

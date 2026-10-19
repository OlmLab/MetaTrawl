# Sync and resume

MetaTrawl's sync commands derive outstanding work from durable project state.
This document explains the checkpoints and incremental behavior used by each
stage.

## Shared principle

Every `sync` command follows the same three rules:

- **Idempotent.** Running it twice with no new inputs is a no-op. It reports what
  is already up to date and exits.
- **Resumable.** If it is interrupted — a crash, a killed job, a cluster
  preemption — rerunning it continues from the last durable checkpoint rather
  than starting over.
- **Incremental.** When you add new SRA runs or local reads, a rerun
  processes only the new work and leaves finished work untouched.

The source of truth for "what is done" is durable state, never in-memory
progress: the DuckDB project store (which runs are `complete`), the per-genome
HDF5 matrix files (which samples they already contain), the comparison DuckDB
files (which sample pairs are already computed), and the shared genome cache
(which genomes are already downloaded and annotated). Scratch space is
disposable; the durable state is authoritative.

## `sync-profile`: profile remaining runs

**Work set.** MetaTrawl profiles the *remaining samples*: every registered input that
is not soft-deleted and does not yet have a `complete` sample. A run marked
`failed` stays in this set, so a rerun automatically retries it. A run that
imported successfully is `complete` and is skipped.

**Per-sample checkpointing.** Each run gets its own scratch directory, and up to
`sample_workers` runs are profiled concurrently. Samples are independent: one
failing sample never blocks the others.

**Per-stage resume.** Within a sample, every stage checks for its own valid
output before running, so an interrupted sample resumes mid-flight:

| Stage | Skips when |
| --- | --- |
| SRA download (`prefetch`, `fasterq-dump`) | non-empty FASTQs already exist in scratch |
| Local-read staging | registered FASTQs have already been copied into scratch |
| Sylph genome selection | an `accessions.txt` already exists |
| Reference preparation (shared cache) | the concatenated per-sample reference is already built |
| Bowtie2 index + alignment | a complete Bowtie2 index and a non-empty BAM already exist |
| ZipStrain `profile-single` | a complete published output bundle already exists |

With fixed-reference profiling, the Sylph and genome-cache rows do not apply.
SRA downloads use `prefetch --max-size u`, so accessions larger than its default
20 GB limit are allowed. MetaTrawl requires a non-empty downloaded archive before
starting `fasterq-dump`; a missing or empty archive fails at the download step
without submitting conversion retries. Rerunning `sync-profile` retries the
download. A retained archive is reused if conversion failed. Ensure enough
scratch space for the archive, FASTQs, and conversion temporary files.

New SRA conversions write to an isolated `.fasterq/` directory under that
sample's scratch, including their temporary files. A pending marker prevents
partially written FASTQs from being reused after interruption. Only successful
conversion publishes the reads and clears the marker; retries overwrite only
the isolated conversion outputs. Existing unmarked read checkpoints remain
usable; this does not retrospectively validate older FASTQs.

Sylph diagnostic fields such as `Eff_lambda` can legitimately contain `HIGH`
or `LOW`. Accession extraction and abundance import preserve those values as
text, convert only required numeric fields, and remove exact repeated headers
with a warning. Other malformed abundance values are errors, not silently
discarded rows. Input tables are never rewritten by these readers.

If ZipStrain exits without a complete, non-empty set of output files, MetaTrawl
does not publish a ready bundle or mark the sample imported. When genome stats
are missing, it reports whether the profile has zero rows or whether the
read-location extraction needs investigation. The row count is read from the
Parquet footer, not by loading the profile. The BAM and workspace are retained.
MetaTrawl never invents genome stats or relaxes read filters to force success;
empty-input/read-location failures in older ZipStrain require diagnosis there.

When SRA conversion produces `sample_1.fastq`, `sample_2.fastq`, and
`sample.fastq`, MetaTrawl aligns the first two as a pair and the third as
unpaired reads, then streams all alignments through one coordinate sort. This
works with Bowtie2 and bwa-mem3, including compressed FASTQs. The sample remains
paired-end: the configured read-inclusion policy is respected, so `proper-pairs`
still excludes unpaired alignments from profiling. Sorting uses a unique
temporary directory for each attempt and cleans it on normal exit or failure.

The prepared ZipStrain assets and Bowtie2 index are shared by every sample and
are reused when their reference/gene/STB content hash is unchanged.

**Commit-then-clean.** When a sample finishes, MetaTrawl imports its ZipStrain
outputs into DuckDB in the coordinator thread (serially, so there is never more
than one DuckDB writer), marks the run `complete`, and only then deletes that
sample's scratch directory and imported output files. If a sample fails, its
scratch is **retained** as a checkpoint and the run is marked `failed` so the
next `sync-profile` retries exactly that run.

**Shared cache.** All workers share one genome/Prodigal cache. A genome that one
sample downloads and annotates is reused by every later sample that needs it,
across runs and across invocations.

## `cache sync-matrix-files`: derive matrix inputs

Matrix building needs a per-genome BED file, STB file, and gene-range table.
MetaTrawl writes these automatically when genomes are downloaded, so most users
never run this command. It exists to (re)derive those files for older caches or
after a forced refresh, straight from the cached `genomes/` and `genes/` FASTAs.
It is idempotent — existing derived files are simply rewritten.

## `matrix sync-build`: one ZipStrain matrix per genome

**Genome set.** By default, every genome represented by at least one `complete`
sample (or just the genomes named with `--genome`).

**Per-genome decision.** For each genome, MetaTrawl looks at `matrices/<genome>.h5`:

- **absent** → build a new ZipStrain matrix from all eligible samples;
- **present** → compare the filtered genome-stat sample set with the compact
  required-sample checkpoint stored in HDF5, then append only missing samples;
- **present, nothing new** → report the genome as up to date.

The matrix keeps committed row order separate from a sorted required-sample
catalog. Count and digest attributes make repeat status checks constant-size;
older matrices are upgraded from their existing sample catalog the first time
they are checked. Under Slurm, MetaTrawl summarizes and checks only one
worker-sized window of genomes at a time, submits the genomes needing work, and
continues checking while those jobs run. It does not materialize every
genome/sample pair or queue every genome before useful work starts.

While an absent matrix is being built, complete sample batches are flushed to a
durable `.h5.tmp` checkpoint. Repeating the same command after cancellation
resumes from the committed sample count; only a batch interrupted during its
write is repeated. The checkpoint becomes the final `.h5` through an atomic
rename when the build completes.

`export_batch_mb` is the approximate RAM target for one matrix batch. MetaTrawl
uses the genome length and matrix dtype to estimate how many complete samples fit,
fetches those samples together from DuckDB, converts them in memory, and appends
the ordered batch to HDF5 with one commit.

**Eligibility.** A sample is eligible for a genome's matrix when it is `complete`
and clears the stat thresholds: `--min-coverage`, `--min-breadth`, and
`--min-ber` (from that genome's ZipStrain genome stats) and `--min-sylph-abundance`
(from the sample's Sylph abundance for that genome). The thresholds are embedded
in the HDF5 file's metadata at build time and **reused automatically on append**,
so every sample added later passes through the same filter as the original build
— you do not repeat the thresholds when appending.

The HDF5 file is the durable handle; the registry is bookkeeping and is not
required for sync behavior.

## `matrix sync-compare`: resumable all-vs-all ANI

For every matrix file in `matrix-dir`, MetaTrawl runs ZipStrain's
`matrix_compare` and writes one comparison DuckDB (`compares/<genome>.duckdb`).
The comparison itself is resumable at the pair level: ZipStrain reopens an
incomplete comparison database and skips pairs that are already computed, so
after `matrix sync-build` appends new samples, a rerun computes only the newly
introduced sample pairs — not the entire matrix. Rerunning with nothing new to
compare is a no-op.

## `sync-genome-views`: static browser bundles

For each completed comparison, MetaTrawl writes a self-contained
`genome_views/<genome>/` bundle (heatmap, dendrogram, clusters, neighbor network,
distributions, and the raw matrices). A bundle is regenerated when either the
comparison or the relevant project statistics change, and skipped otherwise;
older schema-1 bundles are rebuilt as the current schema. Completed pre-1.0
ZipStrain comparison databases are read read-only and never modified.

## Putting it together

Because all five commands key off durable state, the normal way to grow a project
is simply to register more SRA runs and rerun the same five commands. Each one
picks up only its share of the new work:

```bash
metatrawl runs add --db metatrawl.duckdb SRR000010 SRR000011   # add more samples
metatrawl sync-profile        --db metatrawl.duckdb ...          # profiles only the new runs
metatrawl matrix sync-build   --db metatrawl.duckdb ...          # appends only new eligible samples
metatrawl matrix sync-compare --db metatrawl.duckdb ...          # computes only new sample pairs
metatrawl sync-genome-views   --db metatrawl.duckdb ...          # refreshes only changed bundles
```

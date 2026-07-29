# MetaTrawl

<p align="center">
  <img src="metatrawl-concept.svg" alt="MetaTrawl gathers bacterial chromosomes from large metagenomic collections and organizes them for strain-level comparison" width="100%">
</p>

MetaTrawl streamlines [ZipStrain](https://github.com/parsaghadermarzi/ZipStrain)
for very large-scale, strain-level analysis of metagenomic samples. It turns a
collection of SRA run IDs into a durable, queryable project database and
coordinates the expensive steps needed to profile and compare thousands of
samples efficiently.

## MetaTrawl is a wrapper around ZipStrain

MetaTrawl is **not** a replacement for ZipStrain and does not reimplement any of
its science. Every heavy computation in the pipeline is a ZipStrain call that
MetaTrawl prepares inputs for, runs, and files the outputs from. MetaTrawl is the
orchestration and data-management layer that makes running ZipStrain at the scale
of thousands of metagenomes practical.

The division of labor is strict:

| ZipStrain does the science | MetaTrawl does the plumbing |
| --- | --- |
| Profiles a sample against a reference (`profile-single`) | Registers SRA runs, downloads reads, selects the per-sample reference, aligns, and feeds ZipStrain |
| Builds the null model and profiling contract (`prepare_profiling`, `build-null-model`) | Runs those steps once per sample and caches their inputs |
| Computes the strain matrix and all-vs-all ANI (`matrix_pairs.matrix_compare`) | Decides which samples are eligible, assembles the matrix, and stores the comparison |
| Defines the HDF5 matrix and comparison-database formats | Keeps one durable matrix and comparison per genome and appends to them incrementally |

Under the hood, MetaTrawl drives these ZipStrain touchpoints for you:

- `zipstrain utilities prepare_profiling` — turn a reference into a BED file,
  gene-range table, and profiling contract.
- `zipstrain utilities build-null-model` — build the per-sample null model.
- `zipstrain utilities profile-single` — the actual strain profiling step.
- ZipStrain's HDF5 matrix layout — every matrix MetaTrawl writes is a native
  ZipStrain matrix store.
- `zipstrain.matrix_pairs.matrix_compare` — the resumable all-vs-all comparison.

Around those, MetaTrawl orchestrates the external tools ZipStrain expects to
already have inputs from: SRA Toolkit (`prefetch`, `fasterq-dump`), Sylph, NCBI
`datasets`, Prodigal, Bowtie2, and Samtools. You never call any of these by
hand — MetaTrawl does, in the right order, with a shared cache and per-sample
checkpointing.

### The whole pipeline is five commands

Once your project database exists, the entire ZipStrain workflow — from raw SRA
run IDs to browsable per-genome strain comparisons — is five `sync` commands you
can rerun as often as you like:

```bash
metatrawl sync-profile        --db metatrawl.duckdb ...   # download → align → ZipStrain profile → import
metatrawl matrix sync-build   --db metatrawl.duckdb ...   # assemble one ZipStrain matrix per genome
metatrawl matrix sync-compare --db metatrawl.duckdb ...   # ZipStrain all-vs-all ANI per matrix
metatrawl sync-genome-views   --db metatrawl.duckdb ...   # static, browser-ready bundles
metatrawl view genomes        --view-dir genome_views     # explore the results
```

Every `sync` command is idempotent, resumable, and incremental: add more runs at
any time, rerun the same five commands, and each one does only the outstanding
work. The full decision logic behind that behavior is documented in
[How Sync Works (The Logic)](#how-sync-works-the-logic).

## Why MetaTrawl?

Running a strain-level workflow over a few samples is straightforward. Running
the same workflow over thousands of metagenomes introduces a different set of
problems:

- Which SRA runs have already been processed?
- How can workers share downloaded genomes and Prodigal annotations safely?
- How can temporary reads, alignments, and per-sample references be removed
  without losing the durable results?
- How can samples be selected for a genome-specific matrix using coverage,
  breadth, BER, or Sylph abundance?
- How can newly added samples be appended without rebuilding matrices or
  recomputing completed sample pairs?
- How can profile, genome, gene, and abundance data be queried without managing
  thousands of loose files?

MetaTrawl addresses these problems with a mutable DuckDB project store and an
incremental workflow:

1. Register SRA runs.
2. Download and screen reads with Sylph.
3. Reuse a shared genome and Prodigal cache.
4. Align reads and profile samples with ZipStrain.
5. Import profile positions, genome statistics, gene statistics, and Sylph
   abundance into DuckDB.
6. Build genome-specific dense or sparse ZipStrain matrices from eligible
   samples.
7. Run resumable strain-level comparisons and compute only newly introduced
   sample pairs.

Many SRA workers can profile samples concurrently while sharing one prepared
reference cache. Per-sample reads, alignments, concatenated references, and
intermediate outputs live in scratch space and are deleted after successful
import. The DuckDB database, genome cache, matrix stores, and comparison
databases remain durable.

## Installation

### Recommended: Bioconda

For the complete workflow, Conda is the recommended installation method:

```bash
conda install bioconda::metatrawl
```

The Bioconda package installs MetaTrawl together with the external tools used by
the pipeline, including ZipStrain, Sylph, Bowtie2, Samtools, Prodigal, SRA
Toolkit, and NCBI Datasets.

Installing into a dedicated environment keeps the workflow isolated:

```bash
conda create -n metatrawl bioconda::metatrawl
conda activate metatrawl
```

### PyPI

MetaTrawl is also available from PyPI:

```bash
pip install metatrawl
```

The PyPI package installs MetaTrawl and its Python dependencies. It does not
install the external command-line tools required by the complete SRA profiling
workflow. Use Bioconda unless those tools are already installed independently.

### Verify The Installation

Check the installed Python packages and external executables:

```bash
metatrawl test
```

Use the strict check before starting a long workflow. It exits non-zero when a
required dependency is unavailable:

```bash
metatrawl check
```

### Development Installation

```bash
pip install -e ".[test]"
```

This installs the local checkout with the test dependencies. External workflow
tools must still be available in the active environment.

## Tutorial: SRA IDs To Comparisons

This is the normal MetaTrawl workflow. It starts from SRA run IDs and ends with
one comparison database per genome.

### 1. Create An Empty Project

```bash
metatrawl init --db metatrawl.duckdb
```

This creates the DuckDB project store. The database tracks SRA runs, imported
profile rows, genome stats, gene stats, Sylph abundance, and matrix/compare
bookkeeping.

For projects that only need allele-presence comparisons, initialize an
`allele-mask` database:

```bash
metatrawl init \
  --db metatrawl.duckdb \
  --profile-storage allele-mask \
  --profile-min-cov 5
```

ZipStrain profiling is unchanged and still produces its normal full profile.
During import, MetaTrawl stores only the covered-position and allele-set masks
needed by bitmask matrices, plus one shared copy of each cached reference
scaffold. The temporary full profile is deleted only after this transaction
commits.

The threshold is part of the database contract and cannot be changed later.
`allele-mask` supports dense or sparse bitmask matrices, popANI, IBS, and
gene-level popANI. Count matrices, conANI, cosANI, raw A/C/G/T profile queries,
and count reconstruction require normal `full` storage.

### 2. Add SRA Runs

```bash
metatrawl runs add \
  --db metatrawl.duckdb \
  SRR000001 SRR000002 SRR000003
```

Check what is registered:

```bash
metatrawl runs list --db metatrawl.duckdb
```

Check the whole project status:

```bash
metatrawl status --db metatrawl.duckdb
```

The status output is intentionally small: active runs, completed samples,
remaining profiles, profile rows, matrices, and compares.

### 3. Profile Remaining Runs And Import Them

Use `sync-profile` for the high-level profile sync. It finds remaining SRA runs,
downloads reads, runs Sylph, prepares the genome cache, aligns reads, runs
ZipStrain profiling, imports completed outputs into DuckDB, and removes
per-sample scratch/results after successful import.

```bash
metatrawl sync-profile \
  --db metatrawl.duckdb \
  --cache-dir cache \
  --scratch-dir scratch \
  --output-dir outputs \
  --sylph-db /full/path/to/gtdb-r220-c200-dbv1.syldb \
  --threads 16
```

`sync-profile` checkpoints each sample independently. If one sample fails,
successful samples are still imported and cleaned. Failed or incomplete runs stay
pending, so rerunning the same command retries only remaining work.

During profiling, MetaTrawl logs compact lines that work well in terminal and
cluster logs:

```text
METATRAWL sample=SRR123 step=sylph status=done genomes=12 elapsed=4.2s
METATRAWL sample=SRR123 step=cache status=done accessions=10 elapsed=28.9s
METATRAWL sample=SRR123 step=cleanup status=done removed=scratch/SRR123
```

Use an absolute `--sylph-db` path when possible. MetaTrawl validates the file
before launching workers.

### 4. Sync Matrix Requirement Files

When genomes are downloaded, MetaTrawl automatically creates per-genome matrix
requirement files:

```text
cache/genomes/GCF_xxx.fna
cache/genes/GCF_xxx.genes.fna
cache/beds/GCF_xxx.bed
cache/stb/GCF_xxx.stb
cache/gene_ranges/GCF_xxx.gene_ranges.tsv
```

For older caches, or if you want to force-refresh these derived files, run:

```bash
metatrawl cache sync-matrix-files \
  --cache-dir cache
```

For one genome only:

```bash
metatrawl cache sync-matrix-files \
  --cache-dir cache \
  --genome GCF_000269965.1
```

For legacy ZipStrain-style unified reference files, add `--output-dir`:

```bash
metatrawl cache sync-matrix-files \
  --cache-dir cache \
  --output-dir cache/matrix_reference
```

That still writes the per-genome files, plus legacy unified files in
`cache/matrix_reference`.

### 5. Sync Genome Matrices

Build or update one ZipStrain HDF5 matrix per genome represented in the database:

```bash
metatrawl matrix sync-build \
  --db metatrawl.duckdb \
  --matrix-dir matrices \
  --bed-dir cache/beds \
  --stb-dir cache/stb \
  --gene-range-dir cache/gene_ranges \
  --sparse \
  --min-coverage 1 \
  --min-breadth 0.2 \
  --min-ber 0.77 \
  --min-sylph-abundance 0.001
```

For each genome:

- if `matrices/<genome>.h5` does not exist, MetaTrawl builds it;
- if it already exists, MetaTrawl appends eligible samples that are not yet in
  the HDF5 file;
- if no new samples are available, the genome is reported as up to date.

To sync only one genome, add `--genome`:

```bash
metatrawl matrix sync-build \
  --db metatrawl.duckdb \
  --matrix-dir matrices \
  --genome GCF_000269965.1 \
  --bed-dir cache/beds \
  --stb-dir cache/stb \
  --gene-range-dir cache/gene_ranges \
  --sparse
```

The HDF5 matrix file is the durable handle. The old matrix registry is not
required for normal sync behavior.

### 6. Sync Comparisons

Run resumable comparison for every matrix in `matrices/`:

```bash
metatrawl matrix sync-compare \
  --db metatrawl.duckdb \
  --matrix-dir matrices \
  --compare-dir compares \
  --calculate all \
  --backend numpy \
  --memory-limit-gb 16
```

This writes one comparison DuckDB per matrix:

```text
compares/GCF_xxx.duckdb
```

Rerunning `sync-compare` is safe. ZipStrain resumes incomplete comparison
databases and skips completed pairs.

### 7. Sync Genome Views

Prepare self-contained browser-ready artifacts for every completed genome
comparison:

```bash
metatrawl sync-genome-views \
  --db metatrawl.duckdb \
  --compare-dir compares \
  --view-dir genome_views
```

Each `genome_views/<genome>/` bundle is a versioned, self-contained web and
analysis snapshot. A web client starts at `genome_views/catalog.json`, follows
the genome's `manifest.json`, and then loads only the artifacts needed for the
current panel:

- `samples.json`: matrix indices, dendrogram leaf order, missing-data fraction,
  and cluster labels.
- `sample_stats.json`: typed, column-oriented statistics ready for a browser
  table; `sample_stats.parquet` preserves the same data for analytical tools.
- `clusters.json`: reusable clonal and strain assignments, memberships,
  thresholds, and linkage method.
- `dendrogram.json`: the SciPy linkage matrix plus explicit tree merges for an
  interactive dendrogram.
- `neighbor_network.json`: browser-ready nodes and top-neighbor edges, including
  ANI and compared-position counts.
- `similarity_ani.condensed.f32.gz`: the little-endian float32 ANI matrix in
  SciPy-compatible condensed upper-triangle order.
- `total_positions.condensed.u64.gz`: the matching condensed uint64 overlap
  matrix; zero identifies an imputed comparison.
- `distributions.json`: precomputed ANI, overlap, coverage, breadth, BER, and
  abundance histograms.
- `view_data.h5`: the matrices, linkage, ordering, and assignments in a reusable
  scientific container.
- `clustermap.png` and `dendrogram.svg`: static previews, not the primary data
  source for the web interface.

`manifest.json` documents every file's format, media type, size, matrix shape,
dtype, byte order, and axis ordering. The bundle therefore requires no query
against the main MetaTrawl or comparison databases at presentation time.
Rerunning the command skips unchanged bundles, but refreshes them when either
the comparison or relevant project statistics change. Schema-1 bundles are
automatically regenerated as schema 2. Use `--genome GCF_xxx` to refresh one
genome explicitly.

Completed comparison databases from pre-1.0 ZipStrain releases are accepted
read-only. MetaTrawl detects the legacy `genome_pop_ani` result column
automatically. If an older database lacks checkpoint or catalog tables,
MetaTrawl derives completion, samples, and genomes from distinct result rows;
it still skips a view when an available sample catalog shows that result rows
are incomplete. No legacy comparison database is modified.

### Explore Genome Views

Start the interactive genome atlas after `sync-genome-views` completes:

```bash
metatrawl view genomes \
  --view-dir genome_views
```

MetaTrawl serves the generated bundles at `http://127.0.0.1:8766` and opens the
default browser. The viewer provides a searchable genome catalog, summary
distributions, cluster composition, an interactive ANI heatmap, a scalable
dendrogram, a filterable sample-neighbor network, and searchable sample
statistics. It reads only static bundle files and never opens the project or
comparison DuckDB databases.

On a remote cluster, bind to the compute node without attempting to open its
browser:

```bash
metatrawl view genomes \
  --view-dir genome_views \
  --host 0.0.0.0 \
  --port 8766 \
  --no-open
```

Use SSH port forwarding from your workstation:

```bash
ssh -L 8766:127.0.0.1:8766 user@cluster
```

For public deployment, the same viewer and bundles can be hosted as static
files; the local server is not a required production component.

### 8. Inspect Progress

At any point:

```bash
metatrawl status --db metatrawl.duckdb
```

To see which SRA runs still need profile imports:

```bash
metatrawl profiles remaining \
  --db metatrawl.duckdb \
  --output-file remaining_runs.csv
```

## How Sync Works (The Logic)

The tutorial above shows *what* to type. This section explains *why* it is safe
to rerun every command, and how each `sync` decides what work is still
outstanding. Understanding this is the whole point of MetaTrawl: you never track
progress by hand, and you never recompute finished work.

### The shared principle

Every `sync` command follows the same three rules:

- **Idempotent.** Running it twice with no new inputs is a no-op. It reports what
  is already up to date and exits.
- **Resumable.** If it is interrupted — a crash, a killed job, a cluster
  preemption — rerunning it continues from the last durable checkpoint rather
  than starting over.
- **Incremental.** When you add new SRA runs (or new eligible samples), a rerun
  processes only the new work and leaves finished work untouched.

The source of truth for "what is done" is durable state, never in-memory
progress: the DuckDB project store (which runs are `complete`), the per-genome
HDF5 matrix files (which samples they already contain), the comparison DuckDB
files (which sample pairs are already computed), and the shared genome cache
(which genomes are already downloaded and annotated). Scratch space is
disposable; the durable state is authoritative.

### `sync-profile`: profile remaining runs

**Work set.** MetaTrawl profiles the *remaining runs*: every registered run that
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
| Sylph genome selection | an `accessions.txt` already exists |
| Reference preparation (shared cache) | the concatenated per-sample reference is already built |
| Bowtie2 index + alignment | a complete Bowtie2 index and a non-empty BAM already exist |
| ZipStrain `profile-single` | a complete published output bundle already exists |

**Commit-then-clean.** When a sample finishes, MetaTrawl imports its ZipStrain
outputs into DuckDB in the coordinator thread (serially, so there is never more
than one DuckDB writer), marks the run `complete`, and only then deletes that
sample's scratch directory and imported output files. If a sample fails, its
scratch is **retained** as a checkpoint and the run is marked `failed` so the
next `sync-profile` retries exactly that run.

**Shared cache.** All workers share one genome/Prodigal cache. A genome that one
sample downloads and annotates is reused by every later sample that needs it,
across runs and across invocations.

### `cache sync-matrix-files`: derive matrix inputs

Matrix building needs a per-genome BED file, STB file, and gene-range table.
MetaTrawl writes these automatically when genomes are downloaded, so most users
never run this command. It exists to (re)derive those files for older caches or
after a forced refresh, straight from the cached `genomes/` and `genes/` FASTAs.
It is idempotent — existing derived files are simply rewritten.

### `matrix sync-build`: one ZipStrain matrix per genome

**Genome set.** By default, every genome represented by at least one `complete`
sample (or just the genomes named with `--genome`).

**Per-genome decision.** For each genome, MetaTrawl looks at `matrices/<genome>.h5`:

- **absent** → build a new ZipStrain matrix from all eligible samples;
- **present** → append only the eligible samples that are *not already in the
  HDF5 file* (it reads the sample list stored inside the file to diff);
- **present, nothing new** → report the genome as up to date.

**Eligibility.** A sample is eligible for a genome's matrix when it is `complete`
and clears the stat thresholds: `--min-coverage`, `--min-breadth`, and
`--min-ber` (from that genome's ZipStrain genome stats) and `--min-sylph-abundance`
(from the sample's Sylph abundance for that genome). The thresholds are embedded
in the HDF5 file's metadata at build time and **reused automatically on append**,
so every sample added later passes through the same filter as the original build
— you do not repeat the thresholds when appending.

The HDF5 file is the durable handle; the registry is bookkeeping and is not
required for sync behavior.

### `matrix sync-compare`: resumable all-vs-all ANI

For every matrix file in `matrix-dir`, MetaTrawl runs ZipStrain's
`matrix_compare` and writes one comparison DuckDB (`compares/<genome>.duckdb`).
The comparison itself is resumable at the pair level: ZipStrain reopens an
incomplete comparison database and skips pairs that are already computed, so
after `matrix sync-build` appends new samples, a rerun computes only the newly
introduced sample pairs — not the entire matrix. Rerunning with nothing new to
compare is a no-op.

### `sync-genome-views`: static browser bundles

For each completed comparison, MetaTrawl writes a self-contained
`genome_views/<genome>/` bundle (heatmap, dendrogram, clusters, neighbor network,
distributions, and the raw matrices). A bundle is regenerated when either the
comparison or the relevant project statistics change, and skipped otherwise;
older schema-1 bundles are rebuilt as the current schema. Completed pre-1.0
ZipStrain comparison databases are read read-only and never modified.

### Putting it together

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

## Manual Import And Lower-Level Commands

Most users should use `sync-profile`. If an external workflow produced profile
files, import them directly:

```bash
metatrawl profiles import \
  --db metatrawl.duckdb \
  --run-id SRR000001 \
  --profile-file outputs/SRR000001.profile.parquet \
  --genome-stats-file outputs/SRR000001.genome_stats.parquet \
  --gene-stats-file outputs/SRR000001.gene_stats.parquet \
  --sylph-abundance-file outputs/SRR000001.sylph.csv \
  --cache-dir cache
```

Or import many samples from a manifest:

```bash
metatrawl profiles add \
  --db metatrawl.duckdb \
  --manifest completed_profiles.csv \
  --cache-dir cache
```

Manifest columns:

```csv
run_id,profile_file,genome_stats_file,gene_stats_file,sylph_abundance_file
SRR000001,/path/profile.parquet,/path/genome_stats.parquet,/path/gene_stats.parquet,/path/sylph.csv
```

`gene_stats_file` is optional. `--cache-dir` is required only for an
`allele-mask` database when a referenced genome has not already been stored.
A run is complete after its selected profile representation, genome stats, and
Sylph abundance have been imported.

### Compact An Existing Full Database

DuckDB does not return space to the operating system when a large table is
dropped, so conversion writes a new database and never modifies the source:

```bash
metatrawl profiles compact-database \
  --source-db metatrawl.duckdb \
  --output-db metatrawl.allele-mask.duckdb \
  --cache-dir cache \
  --min-cov 5
```

The migration commits one sample at a time. Rerun the same command to resume
after interruption. Do not replace the source until every sample reports
completion and the new database has been validated.

You can still run the lower-level worker command if you want to manage the
remaining-runs CSV yourself:

```bash
metatrawl profiles remaining \
  --db metatrawl.duckdb \
  --output-file remaining_runs.csv

metatrawl profile-sra \
  --db metatrawl.duckdb \
  --remaining-csv remaining_runs.csv \
  --cache-dir cache \
  --scratch-dir scratch \
  --output-dir outputs \
  --sylph-db /full/path/to/gtdb-r220-c200-dbv1.syldb \
  --threads 8
```

## Python Query API

Use a genome view to query one genome across samples, or a sample view to
query everything stored for one sample:

```python
import polars as pl
from metatrawl import open_database

db = open_database("metatrawl.duckdb")

all_genomes = db.genomes().collect()
all_samples = db.samples().collect()
matching_genomes = db.genomes(pattern="996").collect()

genome_stats = db.genome("GCF_000001").genome_stats().collect()
gene_stats = db.genome("GCF_000001").gene_stats().collect()

sample_genomes = db.sample("SRR123").genome_stats().collect()
sample_genes = db.sample("SRR123").gene_stats(genome="GCF_000001").collect()
sample_profile = db.sample("SRR123").profile(genome="GCF_000001").collect()
```

Every query supports `collect()` for a Polars DataFrame, `lazy()` for
additional lazy Polars transformations, and `sink_parquet()` for a direct
DuckDB-to-Parquet export that does not materialize the result in Python:

```python
query = db.genome("GCF_000001").profiles()

query.sink_parquet("GCF_000001.profiles.parquet")
filtered = query.lazy().filter(pl.col("sample_id").is_in(selected_samples))
```

Genome stats, gene stats, Sylph abundance, samples, and genomes remain queryable
in both storage modes. Position-level `profile()` and `profiles()` queries raise
a clear error for `allele-mask` databases because real A/C/G/T counts were
deliberately discarded.

## Controlling concurrency and execution

`sync-profile` and `profile-sra` accept `--workflow-config` with a TOML or JSON file. This separates the number of samples allowed in flight from the concurrency and CPU allocation of each stage. For example, downloads can remain highly parallel while only one Bowtie index build and two alignments run at once:

```bash
metatrawl sync-profile \
  --db metatrawl.duckdb \
  --cache-dir cache \
  --scratch-dir scratch \
  --output-dir outputs \
  --sylph-db /path/to/gtdb.syldb \
  --workflow-config examples/workflow.toml
```

Each stage supports `workers`, `threads`, `execution = "local" | "slurm"`, `retries`, `retry_delay_seconds`, and an optional `environment` table. Slurm stages also accept `time`, `memory_gb`, `partition`, `account`, and arbitrary `extra` `sbatch` options. MetaTrawl submits Slurm jobs with `sbatch --wait`; checkpointing, output publication, and scratch cleanup therefore happen only after the job completes. If a stage fails, MetaTrawl retries that stage command or Slurm job according to the stage retry settings before marking the sample failed.

The configurable stages are `sra_download`, `sylph`, `genome_download`, `prodigal`, `prepare_profile`, `bowtie_build`, `alignment`, `profile`, `matrix_build`, `matrix_compare`, and `genome_view`. Without `--workflow-config`, `--threads` retains the previous profiling behavior.

The same file can configure ZipStrain `profile-single` read filters under `[profile]`: `min_mapq`, `min_baseq`, `min_freq`, `min_read_ani`, and `read_inclusion`. Matrix construction settings belong under `[matrix_build]`: `storage_mode`, `count_dtype`, `min_cov`, `memory_limit_gb`, `export_batch_mb`, and `duckdb_export_threads`. New matrices use compact bitmask storage by default; choose count storage when `conani` or `cosani_<threshold>` is required. An `allele-mask` database rejects count storage and any `min_cov` different from its database contract before creating or resizing a matrix. Comparison settings under `[matrix_compare]` include `calculate`, `ani_method`, `genome`, `backend`, optional `min_cov`, `memory_limit_gb`, and queue/executor controls. When comparison `min_cov` is omitted, ZipStrain uses the build threshold stored in the HDF5 file. `matrix sync-build` can use `[stages.matrix_build]` to submit one build/append job per genome, `matrix sync-compare` can use `[stages.matrix_compare]` to submit one compare job per matrix, and `sync-genome-views` can use `[stages.genome_view]` to submit one artifact job per genome. Explicit CLI values override the matching TOML values.

### Complete TOML template

This copy-ready template configures every pipeline stage. Only alignment uses
Slurm in this example; all other stages run locally. Change a stage's
`execution` to `"slurm"` and give it a corresponding Slurm table when needed.

```toml
# Maximum number of samples progressing through the workflow concurrently.
sample_workers = 12

[stages.sra_download]
workers = 6
threads = 4
execution = "local" # "local" or "slurm"
retries = 2
retry_delay_seconds = 60

[stages.sylph]
workers = 6
threads = 2
execution = "local"
retries = 2
retry_delay_seconds = 60

[stages.genome_download]
workers = 12
threads = 1
execution = "local"
retries = 2
retry_delay_seconds = 60

[stages.prodigal]
workers = 2
threads = 1
execution = "local"
retries = 1
retry_delay_seconds = 30

[stages.prepare_profile]
workers = 4
threads = 2
execution = "local"
retries = 1
retry_delay_seconds = 30

[stages.bowtie_build]
workers = 1
threads = 12
execution = "local"
retries = 1
retry_delay_seconds = 60

[stages.alignment]
workers = 2
threads = 16
execution = "slurm"
retries = 3
retry_delay_seconds = 120

[stages.alignment.slurm]
time = "04:00:00"
memory_gb = 64
partition = "compute"
account = "project-name"

[stages.profile]
workers = 2
threads = 8
execution = "local"
retries = 3
retry_delay_seconds = 120

[stages.matrix_build]
workers = 4
threads = 16
execution = "local"
retries = 1
retry_delay_seconds = 60

[stages.matrix_compare]
workers = 4
threads = 16
execution = "local"
retries = 1
retry_delay_seconds = 60

[stages.genome_view]
workers = 4
threads = 8
execution = "local"
retries = 1
retry_delay_seconds = 60

[profile]
min_mapq = 0
min_baseq = 13
min_freq = 0.0
min_read_ani = 0.95
read_inclusion = "paired"

[matrix_build]
storage_mode = "bitmask" # use "counts" for conANI/cosANI
# count_dtype = "auto"   # valid with storage_mode = "counts"
min_cov = 5
memory_limit_gb = 16
export_batch_mb = 128
duckdb_export_threads = 1

[matrix_compare]
calculate = "all"
ani_method = "popani" # conani and cosani_<threshold> require count matrices
genome = "all"
backend = "numpy"
# min_cov = 5 # normally omit: use the threshold embedded in the matrix
memory_limit_gb = 32
anchor_queue_size = 1
target_queue_size = 2
result_transfer_batch_size = 512
loader_executor_kind = "thread"
writer_executor_kind = "thread"
```

The values above are an example allocation, not universal defaults. Tune
workers, threads, memory, partition, and account for your machine or cluster.
MetaTrawl also accepts optional per-stage `environment` and `slurm.extra` tables
when a real tool or cluster requires them; they are intentionally omitted here
because the standard pipeline does not require any.

Current ZipStrain profiling receives both `reference.fasta` and `profiling_contract.json`. This preserves the reference-aware profile fields and enables `ref_ani` in imported genome and gene statistics. MetaTrawl follows the ZipStrain 1.x defaults of `min_freq = 0`, `min_read_ani = 0.95`, and `read_inclusion = "paired"`; set `min_read_ani = 0` to disable read-ANI filtering. Genuinely single-end inputs continue to use all mapped reads.

MetaTrawl preserves the expanded ZipStrain 1.x statistics in DuckDB. Genome queries include coverage median and standard deviation, genome length, gap statistics, 5x covered sites, heterogeneity, FUG, mapped reads, population and consensus reference ANI, SNS/SNV counts, presence, and taxonomy when those columns are present in the imported ZipStrain output. Gene queries also retain gene length. Older MetaTrawl databases are migrated in place the next time a writable command opens them.

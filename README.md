# MetaTrawl

<p align="center">
  <img src="metatrawl-concept.svg" alt="MetaTrawl gathers bacterial chromosomes from large metagenomic collections and organizes them for strain-level comparison" width="100%">
</p>

MetaTrawl streamlines [ZipStrain](https://github.com/parsaghadermarzi/ZipStrain)
for very large-scale, strain-level analysis of metagenomic samples. It turns a
collection of SRA run IDs into a durable, queryable project database and
coordinates the expensive steps needed to profile and compare thousands of
samples efficiently.

MetaTrawl is not a replacement for ZipStrain. It is the orchestration and data
management layer around it. ZipStrain performs strain profiling and matrix
comparison; MetaTrawl manages the samples, reference cache, parallel profiling,
stored outputs, matrix membership, and incremental analysis.

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

## Install For Development

```bash
pip install -e ".[test]"
```

MetaTrawl checks for the external tools used by the full workflow:
`zipstrain`, `sylph`, `samtools`, `bowtie2`, `prefetch`, `fasterq-dump`,
`datasets`, and `prodigal`.

```bash
metatrawl test
```

Use the strict checker before long jobs. It exits non-zero if anything required
is missing:

```bash
metatrawl check
```

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

### 7. Inspect Progress

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
  --sylph-abundance-file outputs/SRR000001.sylph.csv
```

Or import many samples from a manifest:

```bash
metatrawl profiles add \
  --db metatrawl.duckdb \
  --manifest completed_profiles.csv
```

Manifest columns:

```csv
run_id,profile_file,genome_stats_file,gene_stats_file,sylph_abundance_file
SRR000001,/path/profile.parquet,/path/genome_stats.parquet,/path/gene_stats.parquet,/path/sylph.csv
```

`gene_stats_file` is optional. A run is complete after profile positions, genome
stats, and Sylph abundance have been imported.

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

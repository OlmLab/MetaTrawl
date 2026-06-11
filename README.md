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

## Database Workflow

Initialize a project database:

```bash
metatrawl init --db metatrawl.duckdb
```

Add SRA run IDs:

```bash
metatrawl runs add --db metatrawl.duckdb SRR000001 SRR000002
```

Export only runs that are not yet complete:

```bash
metatrawl profiles remaining \
  --db metatrawl.duckdb \
  --output-file remaining_runs.csv
```

The CSV contains one column:

```csv
run_id
SRR000001
SRR000002
```

Run the high-level sync. This gets remaining runs from DuckDB, runs the SRA
profiling lifecycle, finds completed profile outputs, imports them into DuckDB,
deletes the imported per-sample files, and logs each step:

```bash
metatrawl sync \
  --db metatrawl.duckdb \
  --cache-dir cache \
  --scratch-dir scratch \
  --sylph-db /path/to/gtdb-r220-c200-dbv1.syldb \
  --output-dir outputs \
  --threads 16
```

MetaTrawl runs `sylph profile` for each sample, saves the abundance table as
`SRR000001.sylph.tsv`, extracts nonzero `GCF_...`/`GCA_...` accessions from it,
and asks the shared cache to prepare those genomes. It then runs
`zipstrain utilities prepare_profiling`, builds a Bowtie2 index, aligns reads
with `bowtie2 | samtools`, and runs `zipstrain utilities profile-single` to
produce the profile and stats files that are imported into DuckDB.

Use an absolute `--sylph-db` path when possible. MetaTrawl validates the file
before launching SRA workers.

`--accessions-dir` is still available as a manual override. It should contain
one accession list per run, produced by Sylph or another genome preselection
step:

```text
accessions/SRR000001.accessions.txt
accessions/SRR000002.accessions.csv
```

Each file can be a plain one-accession-per-line text file or a CSV with an
`accession` column.

`sync` expects per-run outputs in `--output-dir` using these conventional names:

```text
SRR000001.profile.parquet
SRR000001.genome_stats.parquet
SRR000001.gene_stats.parquet      # optional
SRR000001.sylph.csv               # csv, tsv, or parquet
```

After a successful import, these per-sample outputs are removed because DuckDB is
the durable project store. The durable cache is left intact. Use
`--keep-profile-outputs` only when debugging a failed or suspicious run.

After profiling, import completed outputs into DuckDB tables:

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

## Cache Workflow

Prepare one sample reference from accessions using a shared cache:

```bash
metatrawl cache prepare \
  --cache-dir cache \
  --accessions accessions.csv \
  --output-dir scratch/SRR000001/reference
```

For parallel workers, start a local cache server:

```bash
metatrawl cache serve \
  --cache-dir cache \
  --host 127.0.0.1 \
  --port 8765
```

The cache keeps only durable per-accession files:

```text
cache/genomes/GCF_xxx.fna
cache/genes/GCF_xxx.genes.fna
```

Per-sample concatenated references are scratch outputs and should be deleted
after import.

## SRA Worker Lifecycle

`profile-sra` wires the worker lifecycle around remaining runs and scratch
cleanup:

```bash
metatrawl profile-sra \
  --db metatrawl.duckdb \
  --remaining-csv remaining_runs.csv \
  --cache-dir cache \
  --scratch-dir scratch \
  --threads 8
```

Long-running steps emit compact cluster-friendly logs:

```text
METATRAWL sample=SRR123 step=sylph status=done genomes=12 elapsed=4.2s
METATRAWL sample=SRR123 step=cache status=done accessions=10 elapsed=28.9s
METATRAWL sample=SRR123 step=cleanup status=done removed=scratch/SRR123
```

## Matrix Workflow

Build reusable matrix inputs from the genome and Prodigal gene cache:

```bash
metatrawl cache build-matrix-files \
  --genome-dir cache/genomes \
  --gene-dir cache/genes \
  --output-dir cache/matrix_reference
```

For one genome only:

```bash
metatrawl cache build-matrix-files \
  --genome-dir cache/genomes \
  --gene-dir cache/genes \
  --genome GCF_000269965.1 \
  --output-dir cache/matrix_reference/GCF_000269965.1
```

Build a ZipStrain matrix from complete DuckDB samples. Thresholds are applied
before temporary profile parquets are exported:

```bash
metatrawl matrix build \
  --db metatrawl.duckdb \
  --genome GCF_000269965.1_ASM26996v1_genomic.fna \
  --bed-file cache/matrix_reference/genomes_bed_file.bed \
  --stb-file cache/matrix_reference/reference.stb \
  --gene-range-table cache/matrix_reference/gene_range_table.tsv \
  --output-file matrices/binfantis.h5 \
  --min-coverage 1 \
  --min-breadth 0.2 \
  --min-ber 0.77 \
  --min-sylph-abundance 0.001 \
  --sparse
```

Append newly imported complete samples to a registered matrix:

```bash
metatrawl matrix append \
  --db metatrawl.duckdb \
  --matrix-id binfantis
```

Compare a registered matrix:

```bash
metatrawl matrix compare \
  --db metatrawl.duckdb \
  --matrix-id binfantis \
  --output-file compares/binfantis.duckdb \
  --calculate all
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

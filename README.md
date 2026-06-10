# MetaTrawl

MetaTrawl is a mutable DuckDB project store for SRA-scale ZipStrain projects.
It tracks run IDs, imports completed ZipStrain/Sylph outputs into real database
tables, coordinates shared genome cache preparation, and builds ZipStrain matrix
stores from selected samples.

The core idea is simple: many SRA workers can run in parallel, but one cache
owner prepares genome and Prodigal outputs for shared accessions. Workers create
per-sample concatenated references in scratch space, profile the sample, import
the final tables into DuckDB, and then delete scratch files.

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
and asks the shared cache to prepare those genomes.

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

Build a ZipStrain matrix from complete DuckDB samples. Thresholds are applied
before temporary profile parquets are exported:

```bash
metatrawl matrix build \
  --db metatrawl.duckdb \
  --genome GCF_000269965.1_ASM26996v1_genomic.fna \
  --bed-file reference/genomes.bed \
  --stb-file reference/genomes.stb \
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

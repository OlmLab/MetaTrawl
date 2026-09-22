# Getting started

This guide takes a MetaTrawl project from registered reads to matrices,
comparisons, and browser-ready genome views. For execution settings and Slurm
examples, see [Workflow configuration](workflow-configuration.md). For the
checkpoint model behind each command, see [Sync and resume](sync-and-resume.md).

## 1. Create an empty project

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

## 2. Register sample inputs

For SRA inputs, provide a CSV with the standard `Run` field. Other
SRA/ENA metadata columns are allowed and ignored by registration:

```csv
Run,study_accession,bioproject
SRR000001,SRP000001,PRJNA000001
SRR000002,SRP000001,PRJNA000001
```

```bash
metatrawl samples add-sra-ids \
  --db metatrawl.duckdb \
  --input-file sra_samples.csv
```

For local reads, use the deterministic columns `sample_id`, `read_1`, and
`read_2`. Leave `read_2` empty for a single-end sample:

```csv
sample_id,read_1,read_2
sample_001,/data/sample_001_R1.fastq.gz,/data/sample_001_R2.fastq.gz
sample_002,/data/sample_002.fastq.gz,
```

```bash
metatrawl samples add-reads \
  --db metatrawl.duckdb \
  --input-file local_reads.csv
```

Registration validates the local files and stores their absolute paths in
DuckDB. It does not copy read data. `sync-profile` copies a sample's reads into
disposable scratch only when that sample starts. The original files are never
deleted. On a cluster, the registered paths must be visible from the worker
nodes.

The older `metatrawl runs add/list/delete` commands remain available. Runs
created by older MetaTrawl versions are interpreted as SRA inputs.

Check the whole project status:

```bash
metatrawl status --db metatrawl.duckdb
```

The status output is intentionally small: active runs, completed samples,
remaining profiles, profile rows, matrices, and compares.

## 3. Profile and import remaining samples

Use `sync-profile` for the high-level profile sync. It finds remaining samples,
downloads SRA reads or stages registered local reads, runs Sylph, prepares
the genome cache, aligns reads, runs ZipStrain profiling, imports completed outputs into DuckDB, and removes
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

### Profile against a fixed reference

When every sample should be mapped against the same existing reference, provide
the reference FASTA, its Prodigal-compatible nucleotide gene FASTA, and the STB
scaffold-to-genome mapping together:

```bash
metatrawl sync-profile \
  --db metatrawl.duckdb \
  --cache-dir cache \
  --scratch-dir scratch \
  --output-dir outputs \
  --reference-genome /references/reference.fna \
  --reference-genome-genes /references/reference.genes.fna \
  --reference-stb /references/reference.stb \
  --threads 16
```

This mode skips Sylph, NCBI genome downloading, and Prodigal. MetaTrawl stages
the three inputs into a content-addressed directory under
`cache/fixed_references/`, prepares the ZipStrain profiling assets, and builds
the Bowtie2 index once. All sample workers reuse those immutable files. The
three options are required together; the reference may contain multiple
scaffolds or genomes as long as the STB maps every scaffold correctly and the
gene FASTA follows the Prodigal header contract.

Fixed-reference samples have no rows in `sylph_abundance`. Profile positions,
genome statistics, and gene statistics are imported normally. A later matrix
build using `--min-sylph-abundance` therefore excludes these samples. The
reference, gene, and STB content hashes are part of the resume checkpoint, so
changing any input prevents reuse of BAM/profile checkpoints generated against
the previous reference.

## 4. Sync matrix requirement files

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

## 5. Sync genome matrices

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

New builds checkpoint complete sample batches into
`matrices/<genome>.h5.tmp`. If a local or Slurm job is cancelled, rerun the same
command: MetaTrawl validates the checkpoint, discards only an unfinished batch,
and continues with the remaining samples. The working file is atomically renamed
to `.h5` after every selected sample has been committed. Do not delete the
`.h5.tmp` file when resuming a build.

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

## 6. Sync comparisons

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

## 7. Sync genome views

Prepare self-contained browser-ready artifacts for every completed genome
comparison:

```bash
metatrawl sync-genome-views \
  --db metatrawl.duckdb \
  --compare-dir compares \
  --view-dir genome_views
```

By default, a sample is retained when no more than 20% of its comparisons for
that genome are missing after applying `--min-comp-len`. This genome-relative
rule scales from small to very large matrices; configure it with
`--max-null-fraction` or `[genome_view].max_null_fraction`. Use
`--max-null-samples` only when an explicit absolute limit is desired. DuckDB
computes connectivity from compact integer sample indices before retained pair
values are transferred to Python.

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

## Explore genome views

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

## 8. Inspect progress

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

## Manual import and lower-level commands

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

### Compact an existing full database

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

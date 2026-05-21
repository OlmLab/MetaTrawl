# MetaTrawl

MetaTrawl is a small registry and bookkeeping tool for large ZipStrain projects.
It keeps track of SRA runs, completed ZipStrain/Sylph profile bundles, and
matrix stores without forcing the heavy profiling work to run inside MetaTrawl.

The v1 workflow is intentionally simple:

1. Add SRA run IDs to a mutable DuckDB registry.
2. Export the remaining unprofiled run IDs as a CSV.
3. Profile those runs with an external SRA/ZipStrain/Sylph workflow.
4. Import completed profile bundles back into MetaTrawl.
5. Build and compare ZipStrain matrix stores from the registered profiles.

## Install for development

```bash
pip install -e ".[test]"
```

## Basic usage

Initialize a registry:

```bash
metatrawl init --db metatrawl.duckdb
```

Add SRA runs:

```bash
metatrawl runs add --db metatrawl.duckdb SRR000001 SRR000002
```

Write the remaining runs that still need profiling:

```bash
metatrawl profiles remaining \
  --db metatrawl.duckdb \
  --output-file remaining_runs.csv
```

The remaining CSV has one column:

```csv
run_id
SRR000001
SRR000002
```

After an external workflow creates ZipStrain and Sylph outputs, import them:

```bash
metatrawl profiles add \
  --db metatrawl.duckdb \
  --manifest completed_profiles.csv
```

The completed profile manifest must contain:

```csv
run_id,profile_file,genome_stats_file,sylph_abundance_file
SRR000001,/path/SRR000001.parquet,/path/SRR000001_genome_stats.parquet,/path/SRR000001_sylph.tsv
```

Build a matrix store for one genome across all complete profiles:

```bash
metatrawl matrix build \
  --db metatrawl.duckdb \
  --genome GCF_000269965.1_ASM26996v1_genomic.fna \
  --output-file matrices/binfantis.h5
```

Compare a registered matrix store:

```bash
metatrawl matrix compare \
  --db metatrawl.duckdb \
  --matrix-id binfantis \
  --output-file compares/binfantis.duckdb \
  --calculate all
```

Check local tool availability:

```bash
metatrawl test
```

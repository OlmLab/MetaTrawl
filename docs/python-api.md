# Python query API

MetaTrawl exposes genome-centered and sample-centered queries over the project
database. Queries can be collected as Polars DataFrames, extended lazily, or
written directly to Parquet.

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

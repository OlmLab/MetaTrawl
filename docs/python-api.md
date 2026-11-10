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

## SNPs and SNVs relative to the reference

`snps()` works in both `full` and `allele-mask` databases. Use it to find
positions where a sample has at least one non-reference allele:

```python
# One sample and genome, as a Polars DataFrame.
snps = db.sample("SRR123").snps("GCF_000001").collect()

# One genome across samples, optionally scoped to a scaffold.
db.genome("GCF_000001").snps(chrom="contig_1").sink_parquet("variants.parquet")

# The same filters can be combined at the database level.
query = db.snps(sample_id="SRR123", genome="GCF_000001")
```

Both modes return one row per sample and variant position, with `sample_id`,
`genome`, `chrom`, one-based scaffold coordinate `pos`, `ref_base`,
`alt_bases` (a list of non-reference bases), `allele_count` (all observed
alleles), and `ref_present`. Reference-matching sites and sites with unknown
reference bases are omitted. These are sample-level calls, not a cohort-wide
polymorphism test. Row order is not guaranteed.

- **Full:** Also returns the adjusted `A`, `C`, `G`, `T` counts,
  `position_coverage`, consensus base `con_base`, `ref_freq`, `con_freq`,
  combined `alt_freq`, and ZipStrain-style `class` (`SNS`, `SNV`, `con_SNV`,
  or `pop_SNV`). The default is `min_cov=5`; you can set a different minimum.
  Full profiles must have retained reference bases to call SNPs.
- **Allele-mask:** Uses the coverage threshold fixed at import. Omit
  `min_cov` or pass the same value. Counts, frequencies, and consensus
  cannot be recovered from allele presence alone.

For large exports, use `sink_parquet()` rather than `collect()` or `lazy()`.
Full mode exports from DuckDB; allele-mask mode decodes compressed blocks
and writes variant rows without collecting the entire result in Python.

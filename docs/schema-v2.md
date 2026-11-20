# Database schema v2

MetaTrawl 2.x requires ZipStrain 2.x and creates schema-v2 project databases.
The package version, project schema version, and ZipStrain profile schema are
separate identities. Existing genome-view files and comparison databases are
not converted by opening a project database.

## Existing projects

- Schema-v1 projects can still be queried with the read-only Python API.
- MetaTrawl 2 refuses writable opens of schema-v1 projects, before changing
  tables or metadata. This guard cannot be bypassed by a provenance override.
- Continue v1 profiling, matrix builds, and comparison resumes with the original
  MetaTrawl/ZipStrain installation. Start a new project database for v2 profiles.
- Unversioned historical registries retain the existing explicit metadata-adoption
  path. Adoption marks old provenance unknown; it does not convert old scientific
  results or make v1 and v2 profiles equivalent.

The generation is part of every new profiling contract. In particular, the
v2 `min_freq` default is `0.01`, and gene statistics and selection follow
ZipStrain 2 semantics. These are not settings to infer retroactively for an
old database. No automatic rewrite of a large historical profile table occurs.

## Stored tables

| Table | Grain and contents |
| --- | --- |
| `profile_positions` | Full mode: sample, genome, scaffold position, adjusted A/C/G/T counts, and optional reference-base mask. No per-position gene column. |
| `allele_mask_profile_blocks` | Allele-mask mode: compressed coverage/presence evidence at the fixed import threshold. |
| `genome_stats` | Sample/genome coverage, breadth, BER, reference ANI, mapped reads, and other supplied summary fields. |
| `gene_stats` | Sample/gene statistics, stable UInt64 `gene_id`, name, coordinates, reference ANI, and six selection fields. |
| `gene_catalog` | One authenticated definition per stable UInt64 gene ID, shared across samples. Includes strand, phase, genetic code, and partial flags. |
| `reference_codons` | One reference triplet per gene and zero-based codon index, shared across samples. |
| `codon_counts` | Optional read-phased evidence: `sample_id`, UInt64 `gene_id`, UInt32 `codon_index`, `codon`, UInt64 `count`. One table for the whole project. |
| `annotation_catalogs` | Imported catalog identities and scientific footer metadata. |
| `sample_profile_metadata` | Sample profile generation and profile/codon scientific metadata, including whether codon evidence was retained. |
| `sylph_abundance` | Sample/genome abundance; may be empty for fixed-reference profiling. |

The existing sample registry, cache records, storage contract, and provenance
tables remain. Gene names are resolved through the shared catalog, never through
a gene column attached to profile positions. Overlapping genes remain independent
catalog entries and retain ZipStrain's overlap-safe statistics.

## Codon retention

With `[profile].collect_codons=false` (the default), MetaTrawl retains the gene
catalog and reference codons, but no observed codon counts; selection fields are
null. With collection enabled, ZipStrain creates read-phased counts and
reference-relative selection. MetaTrawl imports both before cleanup. The choice
does not depend on full versus allele-mask positional storage. It follows
ZipStrain's retention decision and never reconstructs haplotypes from marginals.

The importer authenticates catalog/reference identities, reference-codon content,
gene definitions, read-phased evidence, and the source-profile digest. All
sample evidence, statistics, provenance, and completion status commit in one
transaction. Failed imports roll back; published files remain available to retry.
After commit, per-sample profile, catalog, reference-codon, and observed-codon
Parquets may be removed because their durable contents are in DuckDB.

## Manual imports

Keep the v2 companion files with externally generated profile bundles:

```bash
metatrawl profiles import --db project.duckdb --run-id SRR123 \
  --profile-file SRR123.profile.parquet \
  --genome-stats-file SRR123.genome_stats.parquet \
  --gene-stats-file SRR123.gene_stats.parquet \
  --sylph-abundance-file SRR123.sylph.tsv \
  --gene-info-file SRR123.gene_info.parquet \
  --reference-codons-file SRR123.reference_codons.parquet \
  --codon-counts-file SRR123.codon_counts.parquet
```

Omit the observed-codon file only if it was not collected. Bundle manifests also
accept `gene_info_file`, `reference_codons_file`, and `codon_counts_file` columns.
Keep the profiling provenance sidecar beside the profile for API pair checks.
The catalog and reference-codon files are required for v2 imports even when the
observed-codon table was not collected.

## Analysis scope

Matrix builds continue to use positions or masks. Their `all` comparison means
genome ANI, gene ANI, and max identical stretch, not selection. Pairwise dN/dS
and pN/pS, mosaic regions, and coordinate/significance-resolved identical
stretches are on-demand [Python API](python-api.md) analyses. Gene selection
reads retained codons only, scoped to two samples and one genome.

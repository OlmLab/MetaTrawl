# Workflow configuration

Use a TOML or JSON workflow file to control stage concurrency, local or Slurm
execution, retries, resource escalation, profile filters, matrix construction,
comparison, and genome-view generation. Explicit CLI options override matching
configuration values.

`sync-profile` and `profile-sra` accept `--workflow-config` with a TOML or
JSON file. This separates the number of samples allowed in flight from the
concurrency and CPU allocation of each stage. For example, downloads can remain
highly parallel while only one reference index build and two alignments run at once:

```bash
metatrawl sync-profile \
  --db metatrawl.duckdb \
  --cache-dir cache \
  --scratch-dir scratch \
  --output-dir outputs \
  --sylph-db /path/to/gtdb.syldb \
  --workflow-config examples/workflow.toml
```

Each stage supports `workers`, `threads`, `execution = "local" | "slurm"`,
`retries`, `retry_delay_seconds`, and an optional `environment` table. Slurm
stages also accept `time`, `memory_gb`, `memory_retry_coefficient`,
`time_retry_coefficient`, `partition`, `account`, and arbitrary `extra` `sbatch`
options. Both retry coefficients default to `1.0`. MetaTrawl increases memory
only after an out-of-memory failure and increases time only after a timeout;
preempted jobs retry with unchanged resources. MetaTrawl submits Slurm jobs with
`sbatch --wait`; checkpointing, output publication, and scratch cleanup happen
after the job completes. Failed stages are retried according to their stage
policy before the sample is marked failed.

`sync-profile` publishes each completed sample atomically and sends it to one
dedicated DuckDB writer. Profiling continues while that writer imports bounded
microbatches; when the queue reaches either its sample or byte limit, submission
pauses rather than filling scratch. Pending published bundles are rediscovered on
the next run. MetaTrawl deletes the SRA archive after validated FASTQ creation,
deletes reads and per-sample aligner indexes after BAM validation, and deletes the BAM and
sample reference after publication. Published profile files are deleted only
after their DuckDB transaction commits.

The configurable stages are `sra_download`, `sylph`, `genome_download`,
`prodigal`, `prepare_profile`, `bowtie_build`, `alignment`, `profile`,
`matrix_build`, `matrix_compare`, and `genome_view`. Without
`--workflow-config`, `--threads` controls the profiling workflow.

The same file configures ZipStrain `profile-single` under
`[profile]`: `backend`, `alignment_engine`, `min_mapq`, `min_baseq`, `min_freq`, `min_read_ani`, and
`read_inclusion`. Matrix construction settings belong under `[matrix_build]`:
`storage_mode`, `count_dtype`, `min_cov`, `memory_limit_gb`, `export_batch_mb`,
and `duckdb_export_threads`. New matrices use compact bitmask storage by default;
choose count storage when `conani` or `cosani_<threshold>` is required. An
`allele-mask` database rejects count storage and any `min_cov` different from its
database contract before creating or resizing a matrix.

Comparison settings under `[matrix_compare]` include `calculate`, `ani_method`,
`genome`, `backend`, optional `min_cov`, `memory_limit_gb`, and queue/executor
controls. When comparison `min_cov` is omitted, ZipStrain uses the build
threshold stored in the HDF5 file. `matrix sync-build` can use
`[stages.matrix_build]` to submit one build or append job per genome. Likewise,
`matrix sync-compare` uses `[stages.matrix_compare]` per matrix, and
`sync-genome-views` uses `[stages.genome_view]` per genome. Explicit CLI values
override matching TOML values.

## Complete TOML template

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
memory_retry_coefficient = 1.5
time_retry_coefficient = 1.25
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
backend = "python"
alignment_engine = "bowtie2"
min_mapq = 0
min_baseq = 13
min_freq = 0.0
min_read_ani = 0.95
read_inclusion = "paired"

[profile.null_model]
# Optional overrides; omitted fields use installed ZipStrain defaults.
# error_rate = 0.001
# max_total_reads = 50000
# p_threshold = 0.000001
# model_type = "poisson"

[profile_import]
# One larger-than-limit sample is always allowed through by itself.
queue_max_samples = 8
queue_max_gb = 32
batch_max_samples = 4
batch_max_gb = 8
batch_wait_seconds = 1

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

[genome_view]
min_comp_len = 10000
impute_ani = 97.0
max_null_fraction = 0.20
# max_null_samples = 500 # optional absolute override
linkage_method = "average"
neighbor_k = 20
clonal_cluster_threshold = 99.93
strain_cluster_threshold = 99.8
```

The values above are an example allocation, not universal defaults. Tune
workers, threads, memory, partition, and account for your machine or cluster.
MetaTrawl also accepts optional per-stage `environment` and `slurm.extra` tables
when a real tool or cluster requires them; they are intentionally omitted here
because the standard pipeline does not require any.

Current ZipStrain profiling receives both `reference.fasta` and
`profiling_contract.json`. This preserves reference-aware profile fields and
enables `ref_ani` in imported genome and gene statistics. MetaTrawl follows the
ZipStrain 1.x defaults of `min_freq = 0`, `min_read_ani = 0.95`, and
`read_inclusion = "paired"`; set `min_read_ani = 0` to disable read-ANI
filtering. Single-end inputs use all mapped reads.

MetaTrawl 1.1+ installs the Rust profiler plugin automatically and requires
ZipStrain 1.2.2+. Ensure the same MetaTrawl installation is available in the
environment used by the `profile` stage (including Slurm workers).
Set `[profile].backend = "rust_profiler"`, or pass
`--backend rust_profiler` to `metatrawl sync-profile` or `metatrawl profile-sra`.
The CLI option overrides the TOML/JSON setting; the default is `python`.
Other installed ZipStrain entry-point backend names are also accepted.

Backend selection is recorded in execution provenance, but does not change the
scientific profiling contract. Existing compatible outputs are still reused;
switching backends does not force already completed samples to be reprofiled.
This option is independent of `[matrix_compare].backend` and does not change
stage thread counts or local/Slurm execution policy.

### Null-model settings

Use `[profile.null_model]` in TOML, or the equivalent nested JSON object, to
override ZipStrain's sequencing-error model without changing profiling code:

```toml
[profile.null_model]
error_rate = 0.001
max_total_reads = 50000
p_threshold = 0.000001
model_type = "poisson"
```

These example values match ZipStrain 1.2.2. Omit any field to inherit its value
from the installed ZipStrain version. `error_rate` and `p_threshold` must be
finite and strictly between 0 and 1; `max_total_reads` is a positive integer
coverage ceiling, not a limit on total reads in a sample. Only `poisson` is
currently supported.

Both `sync-profile` and `profile-sra` accept the corresponding CLI options:
`--null-model-error-rate`, `--null-model-max-total-reads`,
`--null-model-p-threshold`, and `--null-model-type`. Each supplied option
overrides only its corresponding TOML/JSON field. The same values reach
per-sample preparation, fallback model generation, and fixed-reference assets.

Effective values are part of the scientific contract and fixed-reference cache
key. Leaving all overrides unset preserves existing default contracts/cache
identities. Changing the model does not silently reuse or relabel old profiles,
models, or BAM checkpoints; follow the database-provenance guards and use new
scratch for changed scientific settings.

### Alignment engines

Bowtie2 remains the default, with its existing sensitive end-to-end settings.
To use **bwa-mem3**, install its executable separately and select:

```toml
[profile]
alignment_engine = "bwa-mem3"
```

Alternatively pass `--alignment-engine bwa-mem3` to `sync-profile` or
`profile-sra`; CLI selection overrides TOML/JSON. MetaTrawl does not install
bwa-mem3. Its executable must be in the `PATH` of both indexing and alignment
stages, including Slurm workers. Missing executables produce an actionable
stage error, not a silent fallback to Bowtie2.

The existing `[stages.bowtie_build]` resource settings control reference indexing
for **either** engine; `[stages.alignment]` controls mapping. Both honor thread
counts, execution environments, retries, and Slurm resource settings. bwa-mem3
uses `index -t N` and `mem -t N` with normal defaults, not `--fast`.
Both engines retain mapped records (`samtools view -F 4`) and produce validated,
coordinate-sorted BAMs before ZipStrain profiling. SAM pipeline failures are
propagated even when sorting itself succeeds.

Aligners can produce different mappings (bwa-mem3 defaults to local alignment).
Unlike a profiling backend, alignment-engine choice is part of the scientific
contract. Changing it is rejected for incompatible database imports and BAM
checkpoints; use a new database/scratch directory for a separate analysis or
explicitly acknowledge mixed database provenance with
`--allow-incompatible-profiles`. That flag never makes an old BAM reusable
with new alignment settings. Existing Bowtie2 contracts and shared reference
caches retain their identities.

Developers can add adapters by subclassing `metatrawl.alignment.AlignmentEngine`
and registering a class under the `metatrawl.alignment_engines` entry-point group.
Adapters define index commands/sidecars, SAM commands, and scientific settings;
MetaTrawl keeps execution and BAM handoff independent of the chosen tool.

MetaTrawl preserves expanded ZipStrain 1.x statistics in DuckDB. Genome
queries include coverage median and standard deviation, genome length, gap
statistics, 5x covered sites, heterogeneity, FUG, mapped reads, population and
consensus reference ANI, SNS/SNV counts, presence, and taxonomy when present in
the imported output. Gene queries also retain gene length. Older MetaTrawl
databases are migrated in place the next time a writable command opens them.

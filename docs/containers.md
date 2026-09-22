# Containers

MetaTrawl publishes separate Linux AMD64 images for portable CPU execution and
NVIDIA CUDA comparison.

| Image | PyTorch | Intended use |
| --- | --- | --- |
| `parsaghadermazi/metatrawl:1.0.0` | CPU | Profiling, matrix building, views, and CPU comparison |
| `parsaghadermazi/metatrawl:latest` | CPU | Latest stable CPU image |
| `parsaghadermazi/metatrawl:1.0.0-cuda` | CUDA 12.8 | The complete workflow with GPU comparison support |

All images include MetaTrawl, ZipStrain, Sylph, Bowtie2, Samtools, Prodigal,
SRA Toolkit, NCBI Datasets, HDF5 support, and their Python dependencies.

## CPU image

Pull and verify the stable image:

```bash
docker pull parsaghadermazi/metatrawl:1.0.0
docker run --rm parsaghadermazi/metatrawl:1.0.0 check
```

MetaTrawl uses `/work` as its container working directory. Mount the project
directory there and use the host user ID so generated files remain writable by
the host account:

```bash
docker run --rm \
  --user "$(id -u):$(id -g)" \
  -v "$PWD:/work" \
  parsaghadermazi/metatrawl:1.0.0 \
  init --db metatrawl.duckdb
```

The image entrypoint is `metatrawl`, so commands begin with the MetaTrawl
subcommand (`init`, `sync-profile`, or `matrix`), not a second `metatrawl`.

## NVIDIA CUDA image

The CUDA image requires a Linux NVIDIA host with a driver that supports CUDA
12.8 and the NVIDIA Container Toolkit configured for Docker. Docker Desktop on
macOS cannot pass an NVIDIA GPU into this Linux container.

Pull the image and verify that PyTorch can see the GPU:

```bash
docker pull parsaghadermazi/metatrawl:1.0.0-cuda

docker run --rm --gpus all \
  --entrypoint python \
  parsaghadermazi/metatrawl:1.0.0-cuda \
  -c 'import torch; assert torch.cuda.is_available(); print(torch.cuda.get_device_name(0))'
```

Run one CUDA comparison while mounting the current project:

```bash
docker run --rm --gpus all \
  --user "$(id -u):$(id -g)" \
  -v "$PWD:/work" \
  parsaghadermazi/metatrawl:1.0.0-cuda \
  matrix compare \
  --db metatrawl.duckdb \
  --matrix-file matrices/GCF_000012825.1.h5 \
  --output-file compares/GCF_000012825.1.duckdb \
  --calculate ani+gene \
  --ani-method popani \
  --backend torch-cuda \
  --memory-limit-gb 10 \
  --no-register
```

`--memory-limit-gb` controls the comparison working set. Leave enough GPU
memory for PyTorch and CUDA allocations outside that limit.

To select one GPU on a multi-GPU host:

```bash
docker run --rm --gpus '"device=1"' \
  --user "$(id -u):$(id -g)" \
  -v "$PWD:/work" \
  parsaghadermazi/metatrawl:1.0.0-cuda \
  matrix compare ... --backend torch-cuda
```

## Apptainer on HPC

Clusters commonly expose containers through Apptainer rather than Docker. Pull
the image once on a login node:

```bash
apptainer pull metatrawl-1.0.0-cuda.sif \
  docker://parsaghadermazi/metatrawl:1.0.0-cuda
```

Use `--nv` inside a GPU allocation and bind the project directory:

```bash
apptainer exec --nv \
  --bind "$PWD:/work" \
  --pwd /work \
  metatrawl-1.0.0-cuda.sif \
  metatrawl matrix compare \
  --db metatrawl.duckdb \
  --matrix-file matrices/GCF_000012825.1.h5 \
  --output-file compares/GCF_000012825.1.duckdb \
  --calculate ani+gene \
  --ani-method popani \
  --backend torch-cuda \
  --memory-limit-gb 10 \
  --no-register
```

Verify GPU access in the same allocation when diagnosing a cluster setup:

```bash
apptainer exec --nv metatrawl-1.0.0-cuda.sif \
  python -c 'import torch; print(torch.cuda.is_available(), torch.cuda.get_device_name(0))'
```

## Build locally

Build the CPU image from the repository root:

```bash
docker build --platform linux/amd64 -t metatrawl:local .
```

Build with the official PyTorch CUDA 12.8 wheel index:

```bash
docker build --platform linux/amd64 \
  --build-arg TORCH_INDEX_URL=https://download.pytorch.org/whl/cu128 \
  -t metatrawl:local-cuda .
```

FROM continuumio/miniconda3:latest AS toolchain

ARG PYTHON_VERSION=3.12
ARG TORCH_INDEX_URL=https://download.pytorch.org/whl/cpu
ARG TORCH_PACKAGE=torch
ARG BWAMEM3_VERSION=0.14.0

ENV PIP_NO_CACHE_DIR=1 \
    PIP_DEFAULT_TIMEOUT=180 \
    PIP_RETRIES=10 \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

RUN printf '%s\n' \
        'channels:' \
        '  - conda-forge' \
        '  - bioconda' \
        'channel_priority: strict' \
    > /root/.condarc \
    && conda install --yes \
        "python=${PYTHON_VERSION}" \
        pip \
        bowtie2 \
        "bwa-mem3=${BWAMEM3_VERSION}" \
        ncbi-datasets-cli \
        prodigal \
        samtools \
        sra-tools \
        sylph \
    && conda clean --all --yes

# CPU is the portable default. Override TORCH_INDEX_URL with an official CUDA
# wheel index, such as cu128, when building an NVIDIA image.
# Keep the base environment consistent with Conda's Python dependency metadata.
RUN python -m pip install --upgrade pip setuptools wheel 'ruamel.yaml<0.19' \
    && python -m pip install \
        --index-url "${TORCH_INDEX_URL}" \
        "${TORCH_PACKAGE}"

FROM toolchain AS runtime

WORKDIR /opt/metatrawl

COPY pyproject.toml README.md LICENSE ./
COPY src ./src

RUN python -m pip install . 'polars>=1.34,<2' \
    && python -m pip check \
    && metatrawl check \
    && bwa-mem3 version \
    && python -c 'from zipstrain.profile_backends import load_backend; assert load_backend("rust_profiler").API_VERSION == 1'

WORKDIR /work

ENTRYPOINT ["metatrawl"]
CMD ["--help"]

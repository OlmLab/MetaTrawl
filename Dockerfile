FROM continuumio/miniconda3:latest

ARG PYTHON_VERSION=3.12
ARG TORCH_INDEX_URL=https://download.pytorch.org/whl/cpu
ARG TORCH_PACKAGE=torch

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
        ncbi-datasets-cli \
        prodigal \
        samtools \
        sra-tools \
        sylph \
    && conda clean --all --yes

# CPU is the portable default. Override TORCH_INDEX_URL with an official CUDA
# wheel index, such as cu128, when building an NVIDIA image.
RUN python -m pip install --upgrade pip setuptools wheel \
    && python -m pip install \
        --index-url "${TORCH_INDEX_URL}" \
        "${TORCH_PACKAGE}"

WORKDIR /opt/metatrawl

COPY pyproject.toml README.md LICENSE ./
COPY src ./src

RUN python -m pip install . \
    && metatrawl check

WORKDIR /work

ENTRYPOINT ["metatrawl"]
CMD ["--help"]

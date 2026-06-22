# Conda packaging for MetaTrawl

The PyPI package only pins MetaTrawl's *Python* dependencies. MetaTrawl also
shells out to external bioinformatics executables that pip cannot install. This
directory holds a single Bioconda recipe (`metatrawl/`) that resolves
**everything** at `conda install` time:

| Layer | Provided by | Comes from |
|-------|-------------|------------|
| Python libs (`click`, `duckdb`, `polars`, `rich`) | `metatrawl` run deps | conda-forge |
| `zipstrain` executable + its libs (incl. `pytorch`, `h5py`) | `zipstrain` | bioconda (already published, 0.10.2) |
| `sylph`, `samtools`, `bowtie2`, `prodigal` | `metatrawl` run deps | bioconda |
| `prefetch`, `fasterq-dump` | `sra-tools` | bioconda |
| `datasets` | `ncbi-datasets-cli` | bioconda |

The required-executable list is taken verbatim from
`src/metatrawl/healthcheck.py` (`REQUIRED_EXECUTABLES` and `REQUIRED_PACKAGES`).

`zipstrain` is already on Bioconda (`bioconda::zipstrain` 0.10.2), so MetaTrawl
just lists it as a run dependency — only **one** recipe to submit.

The recipe is `noarch: python`, built from the PyPI sdist:

- `metatrawl` 0.1.5 — `sha256 3ec0a9d37732d979865b029255a93a214f381d88144864b97055c3c604fba90c`
  (hash of the locally built sdist; valid if you publish that exact artifact)

## Blockers to clear before submitting

1. **Release MetaTrawl as 0.1.5.** The MIT `LICENSE` and the
   `license`/`license-files` fields are now in the repo, but the **0.1.3 sdist
   already on PyPI does not contain `LICENSE`**. Bioconda needs the
   `license_file` inside the sdist, so a fresh release is required (Poetry):
   ```bash
   cd MetaTrawl
   rm -rf dist
   poetry build --format sdist                          # already built once; rebuild is fine
   tar tzf dist/metatrawl-0.1.5.tar.gz | grep LICENSE   # must print a match
   poetry publish                                       # uploads the existing dist/ artifact
   shasum -a 256 dist/metatrawl-0.1.5.tar.gz            # must match meta.yaml sha256
   ```
   Publish the artifact you hashed; if you rebuild before uploading, refresh the
   `sha256` in `metatrawl/meta.yaml`.
2. **Recipe maintainer.** `extra.recipe-maintainers` is set to
   `parsaghadermarzi` — make sure that's your GitHub handle.

## Submitting to Bioconda

```bash
# Fork + clone bioconda-recipes
git clone https://github.com/<you>/bioconda-recipes
cd bioconda-recipes

mkdir -p recipes/metatrawl
cp /path/to/MetaTrawl/conda-recipe/metatrawl/meta.yaml recipes/metatrawl/
git checkout -b add-metatrawl
git add recipes/metatrawl && git commit -m "Add metatrawl"
git push origin add-metatrawl
# open PR; Bioconda CI lints + builds.
```

### Lint / build locally first (optional)

Requires a working conda (the one on this machine is a broken pyenv shim).

```bash
conda install -n base -c conda-forge bioconda-utils
bioconda-utils lint  recipes config.yml --packages metatrawl
bioconda-utils build recipes config.yml --packages metatrawl
```

Or a plain conda-build smoke test outside the bioconda tree:

```bash
conda install -n base -c conda-forge conda-build
conda build conda-recipe/metatrawl -c conda-forge -c bioconda
```

## After it lands on Bioconda

```bash
conda install -c conda-forge -c bioconda metatrawl
metatrawl test     # runs the health check; every dependency should report OK
```

## If you'd rather not wait for Bioconda review

The same `meta.yaml` works as-is with `conda build` against a personal
Anaconda.org channel:

```bash
conda build conda-recipe/metatrawl -c conda-forge -c bioconda
anaconda upload $(conda build conda-recipe/metatrawl --output)
# users: conda install -c <you> -c bioconda -c conda-forge metatrawl
```

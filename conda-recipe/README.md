# Bioconda recipe

MetaTrawl is already available in Bioconda. New releases update the existing
`recipes/metatrawl/meta.yaml` recipe in
[`bioconda/bioconda-recipes`](https://github.com/bioconda/bioconda-recipes).

The recipe installs the Python package from its PyPI source distribution and
adds the external workflow tools that pip cannot provide: Sylph, Samtools,
Bowtie2, Prodigal, SRA Toolkit, and NCBI Datasets.

MetaTrawl 1.1.1 requires ZipStrain 1.2.2+ and `zipstrain-rust-profiler` as a
regular runtime dependency. The native plugin has its own architecture-specific
recipe in `zipstrain-rust-profiler/`; Conda installs it automatically with
MetaTrawl. Submit that recipe alongside the MetaTrawl update so users do not need
a second install command. bwa-mem3 remains an optional user-installed aligner;
the complete MetaTrawl containers include it.

## Release values

| Field | Value |
| --- | --- |
| Version | `1.1.1` |
| Build number | `0` |
| PyPI sdist | `metatrawl-1.1.1.tar.gz` |
| SHA256 | `11c2b13ef29658457c4c3dfb2517aa30d15d375c2fe92ce8fc8660c5f792f80a` |

Verify the release artifact hash against PyPI before submitting the update.
Do not replace it with the hash of a later locally rebuilt archive:

```bash
curl -fsSL https://pypi.org/pypi/metatrawl/<version>/json |
  python -c 'import json,sys; d=json.load(sys.stdin); print(next(x["digests"]["sha256"] for x in d["urls"] if x["packagetype"] == "sdist"))'
```

## Submit the update

Update your fork of `bioconda-recipes`, create a release branch, and replace the
existing recipe:

```bash
git clone https://github.com/<github-user>/bioconda-recipes.git
cd bioconda-recipes
git remote add upstream https://github.com/bioconda/bioconda-recipes.git
git fetch upstream
git checkout master
git rebase upstream/master
git checkout -b update-metatrawl-1.1.1

cp /path/to/MetaTrawl/conda-recipe/metatrawl/meta.yaml \
  recipes/metatrawl/meta.yaml
cp -R /path/to/MetaTrawl/conda-recipe/zipstrain-rust-profiler \
  recipes/zipstrain-rust-profiler

git add recipes/metatrawl/meta.yaml recipes/zipstrain-rust-profiler
git commit -m "Update metatrawl to 1.1.1"
git push -u origin update-metatrawl-1.1.1
```

Open a pull request from that branch to `bioconda/bioconda-recipes:master`.
Because MetaTrawl already has a recipe, do not create a second recipe directory.
Keep `build.number` at `0` for the new upstream version. Increment it only when
rebuilding the same MetaTrawl version with recipe-only corrections.

## Validate locally

From a `bioconda-recipes` checkout with `bioconda-utils` installed:

```bash
bioconda-utils lint recipes config.yml --packages metatrawl zipstrain-rust-profiler
bioconda-utils build recipes config.yml --packages zipstrain-rust-profiler metatrawl
```

After the pull request is merged and the package is published:

```bash
conda create -n metatrawl-1 -c conda-forge -c bioconda metatrawl=1.1.1
conda activate metatrawl-1
metatrawl check
```

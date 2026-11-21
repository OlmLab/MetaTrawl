# Bioconda recipe

MetaTrawl already has a recipe in
[`bioconda/bioconda-recipes`](https://github.com/bioconda/bioconda-recipes/tree/master/recipes/metatrawl).
Update that recipe for each PyPI release rather than adding a second package.

## MetaTrawl 2.0.1

| Field | Value |
| --- | --- |
| Version | `2.0.1` |
| Build number | `0` |
| PyPI sdist | `metatrawl-2.0.1.tar.gz` |
| SHA256 | `10e2fba72ece6954342d9bf53c85f883593aa809546aff89c86e45063140e6b6` |

The runtime requirements must include ZipStrain `>=2.0.0,<3`, its native
Rust profiler `>=0.2.0,<0.3`, and Polars `>=1.36.1,<3`. Preserve these bounds
when updating the recipe; changing only the version and archive hash is not
sufficient for this major release.

The ZipStrain and native-plugin updates must be available first:

- [ZipStrain 2.0.0 update](https://github.com/bioconda/bioconda-recipes/pull/72008)
- [Rust profiler 0.2.0 update](https://github.com/bioconda/bioconda-recipes/pull/70147)

Until they are merged and published, use PyPI or the MetaTrawl 2 containers.
The historical plugin recipe included in this directory is not the v2 recipe;
the plugin's current packaging is maintained with ZipStrain and upstream Bioconda.

## Verify the source artifact

Use the hash published by PyPI, not the hash of a later rebuilt archive:

```bash
curl -fsSL https://pypi.org/pypi/metatrawl/2.0.1/json |
  python -c 'import json,sys; d=json.load(sys.stdin); print(next(x["digests"]["sha256"] for x in d["urls"] if x["packagetype"] == "sdist"))'
```

## Submit and test

On a branch of your Bioconda fork, replace only
`recipes/metatrawl/meta.yaml` with `conda-recipe/metatrawl/meta.yaml`, run the
recipe checks, and open a pull request against `bioconda/bioconda-recipes:master`.
Use build number `0` for a new upstream version. Increment it only for a
recipe-only rebuild of the same version.

```bash
bioconda-utils lint recipes config.yml --packages metatrawl
bioconda-utils build recipes config.yml --packages metatrawl
```

The recipe installs from the source distribution without pip dependency
resolution and declares all runtime dependencies explicitly. Tests cover the
CLI, schema-v2 initialization, pair API imports, and native backend loading.
Bioinformatics tools remain Conda runtime dependencies; bwa-mem3 stays optional.

After the dependencies and MetaTrawl recipe have been published:

```bash
conda create -n metatrawl-2 -c conda-forge -c bioconda metatrawl=2.0.1
conda activate metatrawl-2
metatrawl check
```

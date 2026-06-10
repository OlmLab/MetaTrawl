"""MetaTrawl: registry tooling for SRA-scale ZipStrain profiling projects."""

from importlib.metadata import PackageNotFoundError, version


try:
    __version__ = version("metatrawl")
except PackageNotFoundError:
    __version__ = "0.1.3"

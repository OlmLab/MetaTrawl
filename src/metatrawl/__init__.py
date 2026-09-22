"""MetaTrawl: registry tooling for SRA-scale ZipStrain profiling projects."""

from importlib.metadata import PackageNotFoundError, version

from metatrawl.api import MetaTrawlDatabase, Query, open_database


try:
    __version__ = version("metatrawl")
except PackageNotFoundError:
    __version__ = "1.0.0"


__all__ = ["MetaTrawlDatabase", "Query", "open_database"]

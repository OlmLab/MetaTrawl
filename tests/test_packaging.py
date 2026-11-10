"""Release packaging keeps the native profiler installed but aligners opt-in."""

from pathlib import Path
import tomllib

from packaging.requirements import Requirement


ROOT = Path(__file__).resolve().parents[1]


def test_rust_profiler_is_a_regular_dependency():
    project = tomllib.loads((ROOT / "pyproject.toml").read_text())["project"]
    dependencies = {Requirement(value).name: Requirement(value) for value in project["dependencies"]}
    assert "1.3.0" in dependencies["zipstrain"].specifier
    assert "1.2.3" not in dependencies["zipstrain"].specifier
    assert "0.1.0" in dependencies["zipstrain-rust-profiler"].specifier
    assert "1.34" in dependencies["polars"].specifier
    assert "2.0.0" not in dependencies["polars"].specifier
    assert "bwa-mem3" not in dependencies
    assert "rust" not in project["optional-dependencies"]


def test_container_includes_both_alignment_and_profiling_alternatives():
    dockerfile = (ROOT / "Dockerfile").read_text()
    assert '"bwa-mem3=${BWAMEM3_VERSION}"' in dockerfile
    assert "ARG BWAMEM3_VERSION=0.14.0" in dockerfile
    assert "bwa-mem3 version" in dockerfile
    assert 'load_backend("rust_profiler")' in dockerfile
    assert 'ARG ZIPSTRAIN_PACKAGE="zipstrain>=1.3.0,<2"' in dockerfile
    assert '"${ZIPSTRAIN_PACKAGE}"' in dockerfile
    assert '{"anchor", "block"} <= set(mp.MATRIX_COMPARE_KERNELS)' in dockerfile

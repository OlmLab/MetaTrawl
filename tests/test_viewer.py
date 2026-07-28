from __future__ import annotations

import json
from pathlib import Path
from click.testing import CliRunner
import pytest

from metatrawl import cli
from metatrawl import viewer


def _view_dir(tmp_path: Path) -> Path:
    view_dir = tmp_path / "genome_views"
    view_dir.mkdir()
    (view_dir / "catalog.json").write_text(
        json.dumps(
            {
                "schema_version": 2,
                "genome_count": 1,
                "genomes": [{"genome": "GCF_1", "path": "GCF_1"}],
            }
        )
    )
    bundle = view_dir / "GCF_1"
    bundle.mkdir()
    (bundle / "manifest.json").write_text('{"schema_version": 2, "genome": "GCF_1"}')
    return view_dir


def test_genome_view_handler_routes_application_and_bundle_data_safely(tmp_path: Path) -> None:
    view_dir = _view_dir(tmp_path)
    handler_type = viewer.genome_view_handler(view_dir)
    handler = handler_type.__new__(handler_type)

    assert Path(handler.translate_path("/")).read_text().find("MetaTrawl Genome Atlas") >= 0
    assert Path(handler.translate_path("/app.js")).read_text().find('const DATA = "/data/"') >= 0
    styles = Path(handler.translate_path("/styles.css")).read_text()
    assert '"DINish"' in styles
    assert '"Azeret Mono"' in styles
    assert Path(handler.translate_path("/fonts/DINish-Regular.woff2")).stat().st_size > 10_000
    assert Path(handler.translate_path("/fonts/AzeretMono-Latin.woff2")).stat().st_size > 10_000
    assert json.loads(Path(handler.translate_path("/data/catalog.json")).read_text())["genome_count"] == 1
    assert json.loads(Path(handler.translate_path("/data/GCF_1/manifest.json")).read_text())["genome"] == "GCF_1"
    assert not Path(handler.translate_path("/data/../../pyproject.toml")).exists()


def test_genome_view_directory_requires_catalog(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError, match="sync-genome-views"):
        viewer.validate_genome_view_directory(tmp_path)


def test_view_genomes_command_is_registered() -> None:
    result = CliRunner().invoke(cli.cli, ["view", "genomes", "--help"])
    assert result.exit_code == 0
    assert "--view-dir" in result.output
    assert "--no-open" in result.output

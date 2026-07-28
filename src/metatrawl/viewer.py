"""Static browser application and local server for genome view bundles."""

from __future__ import annotations

from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
import json
from pathlib import Path
from typing import Type
from urllib.parse import unquote, urlsplit


ASSET_DIR = Path(__file__).with_name("web") / "genomes"


def validate_genome_view_directory(view_dir: Path) -> Path:
    """Return a resolved view directory containing a valid catalog."""
    resolved = Path(view_dir).expanduser().resolve()
    catalog = resolved / "catalog.json"
    if not catalog.is_file():
        raise FileNotFoundError(
            f"Genome view catalog does not exist: {catalog}. "
            "Run `metatrawl sync-genome-views` first."
        )
    try:
        payload = json.loads(catalog.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError(f"Genome view catalog is not valid JSON: {catalog}") from exc
    if not isinstance(payload.get("genomes"), list):
        raise ValueError(f"Genome view catalog is missing its genomes list: {catalog}")
    if not (ASSET_DIR / "index.html").is_file():
        raise RuntimeError(f"MetaTrawl genome viewer assets are missing: {ASSET_DIR}")
    return resolved


def genome_view_handler(view_dir: Path) -> Type[SimpleHTTPRequestHandler]:
    """Create an HTTP handler serving application assets and `/data/` bundles."""
    data_root = validate_genome_view_directory(view_dir)
    asset_root = ASSET_DIR.resolve()

    class GenomeViewHandler(SimpleHTTPRequestHandler):
        def translate_path(self, path: str) -> str:
            request_path = unquote(urlsplit(path).path)
            if request_path.startswith("/data/"):
                root = data_root
                relative = request_path.removeprefix("/data/")
            else:
                root = asset_root
                relative = request_path.lstrip("/") or "index.html"
            candidate = (root / relative).resolve()
            try:
                candidate.relative_to(root)
            except ValueError:
                return str(root / "__not_found__")
            return str(candidate)

        def end_headers(self) -> None:
            self.send_header("Cache-Control", "no-cache")
            self.send_header("X-Content-Type-Options", "nosniff")
            self.send_header("Referrer-Policy", "no-referrer")
            super().end_headers()

        def log_message(self, format: str, *args: object) -> None:
            return

    return GenomeViewHandler


def create_genome_view_server(
    view_dir: Path,
    *,
    host: str = "127.0.0.1",
    port: int = 8766,
) -> ThreadingHTTPServer:
    """Create, but do not start, the threaded genome viewer server."""
    return ThreadingHTTPServer((host, port), genome_view_handler(view_dir))

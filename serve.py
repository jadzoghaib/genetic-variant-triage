"""Static server that gzips — a closer match to how the site is actually served.

`python -m http.server` sends no `Content-Encoding`, so BRCA1's variants arrive
as 1.1 MB locally while GitHub Pages would send roughly a quarter of that. That
gap makes local testing pessimistic in a way that hides how the site really
performs, and encourages optimising the wrong thing.

    uv run python serve.py           # http://localhost:8080
    uv run python serve.py --port 9000 --dir site
"""

from __future__ import annotations

import argparse
import gzip
import io
from functools import partial
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

# Only text-ish payloads benefit; .pdb and .json compress extremely well,
# images and already-compressed formats do not.
COMPRESSIBLE = {".json", ".js", ".css", ".html", ".pdb", ".svg", ".map", ".txt"}
MIN_BYTES = 1024


class GzipHandler(SimpleHTTPRequestHandler):
    protocol_version = "HTTP/1.1"

    def end_headers(self):
        # A console is only as trustworthy as the freshness of what it shows;
        # during development a cached payload silently contradicts a rebuild.
        self.send_header("Cache-Control", "no-cache")
        super().end_headers()

    def send_head(self):
        path = Path(self.translate_path(self.path))
        accepts_gzip = "gzip" in self.headers.get("Accept-Encoding", "")
        if not (accepts_gzip and path.is_file()
                and path.suffix.lower() in COMPRESSIBLE
                and path.stat().st_size >= MIN_BYTES):
            return super().send_head()

        raw = path.read_bytes()
        body = gzip.compress(raw, compresslevel=6)
        self.send_response(200)
        self.send_header("Content-Type", self.guess_type(str(path)))
        self.send_header("Content-Encoding", "gzip")
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Vary", "Accept-Encoding")
        # Without a validator the browser cannot revalidate a cached copy, and
        # a rebuilt payload can go on being served from cache — which looks
        # exactly like a code bug.
        self.send_header("Last-Modified", self.date_time_string(
            int(path.stat().st_mtime)))
        self.end_headers()
        return io.BytesIO(body)

    def log_message(self, fmt, *args):
        if "200" not in (args[1] if len(args) > 1 else ""):
            super().log_message(fmt, *args)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--port", type=int, default=8080)
    ap.add_argument("--dir", default="site")
    args = ap.parse_args()

    root = Path(__file__).parent / args.dir
    handler = partial(GzipHandler, directory=str(root))
    print(f"serving {root} with gzip on http://localhost:{args.port}")
    ThreadingHTTPServer(("", args.port), handler).serve_forever()


if __name__ == "__main__":
    main()

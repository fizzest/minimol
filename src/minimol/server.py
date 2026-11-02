"""A small local web server: serves the editor UI and a JSON API.

Only the Python standard library is used. The server is meant to be bound to
127.0.0.1 and protects itself against other web pages in the same browser:

* every API call must be a POST carrying the header ``X-Minimol: 1``. Browsers
  will not send a custom header cross-origin without a CORS preflight, which
  this server never approves;
* the ``Host`` header must name this server (blocks DNS-rebinding attacks);
* only files listed at start-up from the package's ``static`` folder are served.
"""

from __future__ import annotations

import json
import mimetypes
import os
import threading
import traceback
from http import HTTPStatus
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

from rdkit import rdBase

from . import __version__, chem, edit, gxtb

STATIC_DIR = Path(__file__).parent / "static"
MAX_BODY = 64 * 1024 * 1024  # bytes
CSP = (
    "default-src 'self'; script-src 'self'; style-src 'self'; img-src 'self' data: blob:; "
    "connect-src 'self'; object-src 'none'; base-uri 'none'; frame-ancestors 'none'"
)
_TYPES = {
    ".html": "text/html; charset=utf-8",
    ".js": "text/javascript; charset=utf-8",
    ".css": "text/css; charset=utf-8",
    ".svg": "image/svg+xml",
    ".txt": "text/plain; charset=utf-8",
}


def _static_files() -> dict[str, Path]:
    files = {}
    for path in STATIC_DIR.rglob("*"):
        if path.is_file():
            files["/static/" + path.relative_to(STATIC_DIR).as_posix()] = path
    files["/"] = STATIC_DIR / "index.html"
    return files


class MinimolServer(ThreadingHTTPServer):
    daemon_threads = True

    def __init__(self, address, initial_file: str | None = None):
        super().__init__(address, Handler)
        self.static_files = _static_files()
        self.initial_file = initial_file
        # One optimization at a time: g-xTB is multi-threaded on its own.
        self.optimize_lock = threading.Lock()

    @property
    def allowed_hosts(self) -> set[str]:
        port = self.server_address[1]
        return {f"127.0.0.1:{port}", f"localhost:{port}", f"[::1]:{port}"}

    @property
    def url(self) -> str:
        return f"http://127.0.0.1:{self.server_address[1]}/"


class Handler(BaseHTTPRequestHandler):
    server: MinimolServer
    server_version = f"minimol/{__version__}"
    protocol_version = "HTTP/1.1"

    def log_message(self, format, *args):  # noqa: A002 - signature from base class
        if os.environ.get("MINIMOL_DEBUG"):
            super().log_message(format, *args)

    # ------------------------------------------------------------- helpers

    def _send(self, status: int, body: bytes, content_type: str, extra: dict | None = None):
        self.send_response(status)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        self.send_header("X-Content-Type-Options", "nosniff")
        self.send_header("Content-Security-Policy", CSP)
        self.send_header("Referrer-Policy", "no-referrer")
        for key, value in (extra or {}).items():
            self.send_header(key, value)
        self.end_headers()
        if self.command != "HEAD":
            self.wfile.write(body)

    def _json(self, status: int, payload: dict):
        self._send(status, json.dumps(payload).encode(), "application/json")

    def _host_ok(self) -> bool:
        return self.headers.get("Host", "") in self.server.allowed_hosts

    # -------------------------------------------------------------- routes

    def do_GET(self):  # noqa: N802 - name from base class
        if not self._host_ok():
            return self._send(HTTPStatus.FORBIDDEN, b"Forbidden host", "text/plain")
        path = self.path.split("?", 1)[0]
        file = self.server.static_files.get(path)
        if file is None:
            return self._send(HTTPStatus.NOT_FOUND, b"Not found", "text/plain")
        ctype = _TYPES.get(file.suffix) or mimetypes.guess_type(file.name)[0] or "application/octet-stream"
        self._send(HTTPStatus.OK, file.read_bytes(), ctype)

    do_HEAD = do_GET

    def do_POST(self):  # noqa: N802 - name from base class
        if not self._host_ok() or self.headers.get("X-Minimol") != "1":
            self.close_connection = True  # the request body is left unread
            return self._json(HTTPStatus.FORBIDDEN, {"ok": False, "error": "Forbidden"})
        try:
            length = int(self.headers.get("Content-Length", "0"))
        except ValueError:
            length = -1
        if not 0 <= length <= MAX_BODY:
            self.close_connection = True
            return self._json(HTTPStatus.REQUEST_ENTITY_TOO_LARGE, {"ok": False, "error": "Request too large"})
        try:
            body = json.loads(self.rfile.read(length) or b"{}")
            if not isinstance(body, dict):
                raise ValueError
        except ValueError:
            return self._json(HTTPStatus.BAD_REQUEST, {"ok": False, "error": "Request body must be a JSON object"})

        route = ROUTES.get(self.path)
        if route is None:
            return self._json(HTTPStatus.NOT_FOUND, {"ok": False, "error": "Unknown API call"})
        try:
            result = route(self.server, body)
        except chem.ChemError as exc:
            return self._json(HTTPStatus.BAD_REQUEST, {"ok": False, "error": str(exc)})
        except Exception as exc:  # report, keep serving
            traceback.print_exc()
            return self._json(
                HTTPStatus.INTERNAL_SERVER_ERROR,
                {"ok": False, "error": f"Internal error ({type(exc).__name__}): {exc}"},
            )
        self._json(HTTPStatus.OK, {"ok": True, **result})


# ------------------------------------------------------------------ API calls


def api_info(server: MinimolServer, body: dict) -> dict:
    binary = gxtb.find_binary(refresh=bool(body.get("refresh")))
    methods = list(chem.FORCE_FIELDS) + (["g-xTB"] if binary else [])
    return {
        "version": __version__,
        "rdkit": rdBase.rdkitVersion,
        "methods": methods,
        "gxtb": binary,
        "export_formats": list(chem.EXPORT_FORMATS),
        "fragments": edit.fragment_list(),
    }


def api_import(server: MinimolServer, body: dict) -> dict:
    molecule, notes = chem.import_text(body.get("name"), body.get("text"))
    return {"molecule": molecule, "notes": notes}


def api_export(server: MinimolServer, body: dict) -> dict:
    title = body.get("title") if isinstance(body.get("title"), str) else "minimol"
    title = " ".join(title.split())[:70] or "minimol"
    return {"text": chem.export_text(body.get("molecule"), body.get("format"), title)}


def api_smiles(server: MinimolServer, body: dict) -> dict:
    molecule, notes = chem.from_smiles(body.get("smiles"))
    return {"molecule": molecule, "notes": notes}


def api_addh(server: MinimolServer, body: dict) -> dict:
    molecule, notes = chem.add_hydrogens(body.get("molecule"))
    return {"molecule": molecule, "notes": notes}


def api_edit(server: MinimolServer, body: dict) -> dict:
    return edit.apply_edit(body.get("molecule"), body.get("op"), body.get("keep_h", True), body.get("view"))


def api_optimize(server: MinimolServer, body: dict) -> dict:
    method = body.get("method", "MMFF94")
    fixed = body.get("fixed", [])
    with server.optimize_lock:
        if method == "g-xTB":
            return gxtb.optimize(
                body.get("molecule"),
                charge=body.get("charge"),
                uhf=body.get("uhf"),
                fixed=fixed,
                level=body.get("level", "normal"),
                max_cycles=body.get("max_iters"),
            )
        return chem.optimize(body.get("molecule"), method, body.get("max_iters") or 2000, fixed)


def api_initial(server: MinimolServer, body: dict) -> dict:
    path = server.initial_file
    if not path:
        return {"molecule": None}
    try:
        text = Path(path).read_text(errors="replace")
    except OSError as exc:
        raise chem.ChemError(f"Could not open {path}: {exc.strerror}") from None
    molecule, notes = chem.import_text(os.path.basename(path), text)
    return {"molecule": molecule, "notes": notes, "name": os.path.basename(path)}


ROUTES = {
    "/api/info": api_info,
    "/api/import": api_import,
    "/api/export": api_export,
    "/api/smiles": api_smiles,
    "/api/addh": api_addh,
    "/api/edit": api_edit,
    "/api/optimize": api_optimize,
    "/api/initial": api_initial,
}


def make_server(host: str = "127.0.0.1", port: int = 8765, initial_file: str | None = None) -> MinimolServer:
    """Create the server; falls back to a free port if ``port`` is taken."""
    try:
        return MinimolServer((host, port), initial_file)
    except OSError:
        if port == 0:
            raise
        return MinimolServer((host, 0), initial_file)

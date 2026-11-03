"""``minimol --self-test REPORT.json``: check that an installation works.

Used by the Windows build to test the packaged program (before and after
installing it) and handy for support: it checks the chemistry (RDKit), the
editor's files and server, and optionally (``--self-test-window``) opens the
real window and waits for the editor to load in it. Results go to a JSON file
because the Windows program has no console.
"""

from __future__ import annotations

import http.client
import json
import platform
import sys
import threading
import time
import traceback

from . import __version__


def _check(results: dict, name: str, fn) -> None:
    try:
        results[name] = {"ok": True, "detail": fn()}
    except Exception as exc:  # recorded, not raised: the report lists every failure
        results[name] = {"ok": False, "detail": f"{type(exc).__name__}: {exc}", "trace": traceback.format_exc()}


def _chemistry() -> dict:
    from rdkit import rdBase

    from . import chem

    mol, _ = chem.from_smiles("CCO")
    result = chem.optimize(mol, "MMFF94")
    return {"rdkit": rdBase.rdkitVersion, "ethanol_mmff_kcal": round(result["energy"], 4), "converged": result["converged"]}


def _editing() -> dict:
    from . import layout
    from .edit import apply_edit

    methane = apply_edit({"atoms": [], "bonds": []}, {"type": "add_atom", "element": "C", "pos": [0, 0, 0]})["molecule"]
    toluene = apply_edit(methane, {"type": "fragment", "atom": 0, "key": "3"})["molecule"]
    drawing = layout.smiles_2d("c1ccccc1O")["molecule"]
    return {"benzene_atoms": len(toluene["atoms"]), "phenol_2d_atoms": len(drawing["atoms"])}


def _server() -> dict:
    from .server import make_server

    srv = make_server("127.0.0.1", 0)
    threading.Thread(target=srv.serve_forever, daemon=True).start()
    port = srv.server_address[1]
    sizes = {}
    try:
        for path in ("/", "/static/js/app.js", "/static/js/viewer2d.js", "/static/vendor/three.module.js", "/static/vendor/three.core.js", "/static/style.css"):
            conn = http.client.HTTPConnection("127.0.0.1", port, timeout=10)
            conn.request("GET", path, headers={"Host": f"127.0.0.1:{port}"})
            res = conn.getresponse()
            body = res.read()
            conn.close()
            if res.status != 200 or not body:
                raise RuntimeError(f"GET {path} -> HTTP {res.status}")
            sizes[path] = len(body)
        conn = http.client.HTTPConnection("127.0.0.1", port, timeout=10)
        conn.request("POST", "/api/info", body=b"{}", headers={"Host": f"127.0.0.1:{port}", "X-Minimol": "1", "Content-Type": "application/json"})
        info = json.loads(conn.getresponse().read())
        conn.close()
        return {"files": sizes, "methods": info["methods"], "gxtb": info["gxtb"], "gxtb_installable": info["gxtb_installable"]}
    finally:
        srv.shutdown()
        srv.server_close()


def _qt_webengine() -> str:
    from . import desktop

    ok, reason = desktop.available()
    if not ok:
        raise RuntimeError(reason)
    return "available"


def _window(software_gl: bool) -> dict:
    """Open the real desktop window, wait until the editor has started in it, close it."""
    from . import desktop
    from .server import make_server

    desktop._prepare_environment(software_gl)
    from PySide6.QtCore import QEvent
    from PySide6.QtWidgets import QApplication

    app = QApplication.instance() or QApplication(sys.argv[:1])
    srv = make_server("127.0.0.1", 0)
    threading.Thread(target=srv.serve_forever, daemon=True).start()
    win = desktop.create_window(srv.url)
    win.show()
    answer = {}
    deadline = time.time() + 90
    try:
        while time.time() < deadline and "state" not in answer:
            win.js(desktop.PROBE_JS, lambda v: answer.__setitem__("last", v))
            for _ in range(20):
                app.processEvents()
                time.sleep(0.01)
            if answer.get("last") in ("ready", "no-webgl"):
                answer["state"] = answer["last"]
        state = answer.get("state", "timed out")
        if state not in ("ready", "no-webgl"):
            raise RuntimeError(f"the editor did not start in the window ({state})")
        return {"page": "loaded", "webgl": state == "ready", "qt_platform": app.platformName()}
    finally:
        win._close_confirmed = True
        win.close()
        win.release_web()
        QApplication.sendPostedEvents(None, QEvent.Type.DeferredDelete)
        srv.shutdown()
        srv.server_close()


def run(report_path: str, window: bool = False, software_gl: bool = False) -> int:
    results: dict = {
        "minimol": __version__,
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "frozen": bool(getattr(sys, "frozen", False)),
    }
    checks: dict = {}
    _check(checks, "chemistry", _chemistry)
    _check(checks, "editing", _editing)
    _check(checks, "server", _server)
    _check(checks, "qt_webengine", _qt_webengine)
    if window:
        _check(checks, "window", lambda: _window(software_gl))
    results["checks"] = checks
    results["ok"] = all(c["ok"] for c in checks.values())
    with open(report_path, "w", encoding="utf-8") as fh:
        json.dump(results, fh, indent=2)
    return 0 if results["ok"] else 1

"""Smoke tests for the desktop window (Qt WebEngine), run with Qt's offscreen platform."""

import json
import os
import threading
import time

import pytest

pytest.importorskip("PySide6.QtWebEngineWidgets", reason="PySide6 with Qt WebEngine not available", exc_type=ImportError)
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PySide6.QtCore import QEvent, Qt  # noqa: E402
from PySide6.QtTest import QTest  # noqa: E402
from PySide6.QtWidgets import QApplication, QFileDialog, QMessageBox  # noqa: E402

from minimol import clipboard, desktop  # noqa: E402
from minimol.server import make_server  # noqa: E402


@pytest.fixture(scope="module")
def qapp():
    # The test machine may have no GPU; allow WebGL on a software renderer.
    desktop._prepare_environment(software_gl=True)
    app = QApplication.instance() or QApplication([])
    yield app
    # Qt's "offscreen" test platform crashes at exit if the clipboard still
    # holds data (verified: the real X11 platform does not), so empty it here.
    # The app itself never clears the clipboard: what you copy stays available.
    app.clipboard().clear()


@pytest.fixture()
def window(qapp, tmp_path):
    water = tmp_path / "water.xyz"
    water.write_text("3\nwater\nO 0 0 0\nH 0.96 0 0\nH -0.24 0.93 0\n")
    server = make_server("127.0.0.1", 0, str(water))
    server.clipboard = desktop.make_clipboard_bridge()
    threading.Thread(target=server.serve_forever, daemon=True).start()
    win = desktop.create_window(server.url)
    win.show()
    wait_for(qapp, lambda: (js(qapp, win, "document.querySelector('#method') ? document.querySelector('#method').options.length : 0") or 0) > 0)
    yield win
    win._close_confirmed = True
    win.close()
    win.deleteLater()
    QApplication.sendPostedEvents(None, QEvent.Type.DeferredDelete)
    server.shutdown()
    server.server_close()


def wait_for(app, predicate, timeout=60):
    end = time.time() + timeout
    while time.time() < end:
        app.processEvents()
        if predicate():
            return
        time.sleep(0.02)
    raise AssertionError("timed out")


def js(app, win, code):
    box = {}
    win.js(code, lambda value: box.setdefault("v", value))
    end = time.time() + 30
    while "v" not in box and time.time() < end:
        app.processEvents()
        time.sleep(0.01)
    return box.get("v")


def formula(app, win):
    return js(app, win, "window.__minimol ? window.__minimol.M.formula(window.__minimol.model()) : ''")


def test_two_windows_can_run_at_once(qapp, window, tmp_path):
    # Each window has its own in-memory profile, so a second one loads fine.
    server = make_server("127.0.0.1", 0, None)
    threading.Thread(target=server.serve_forever, daemon=True).start()
    second = desktop.create_window(server.url)
    second.show()
    try:
        wait_for(qapp, lambda: (js(qapp, second, "window.__minimol ? 1 : 0") or 0) == 1, timeout=30)
    finally:
        second._close_confirmed = True
        second.close()
        second.deleteLater()
        QApplication.sendPostedEvents(None, QEvent.Type.DeferredDelete)
        server.shutdown()
        server.server_close()


def test_window_loads_editor_and_initial_file(qapp, window):
    wait_for(qapp, lambda: formula(qapp, window) == "H2O")
    assert js(qapp, window, "window.minimolDesktop === true")
    wait_for(qapp, lambda: "water" in window.windowTitle())


def test_ctrl_r_optimizes(qapp, window, monkeypatch):
    wait_for(qapp, lambda: formula(qapp, window) == "H2O")
    js(qapp, window, "document.querySelector('#method').value = 'MMFF94'")
    window.view.setFocus()
    target = window.view.focusProxy() or window.view
    QTest.keyClick(target, Qt.Key.Key_R, Qt.KeyboardModifier.ControlModifier)
    wait_for(qapp, lambda: "kcal/mol" in (js(qapp, window, "document.querySelector('#energy').textContent") or ""))
    # Still the same page (no reload): the undo history has the optimization in it.
    assert js(qapp, window, "window.__minimol.history.stack.length") == 3


def test_save_dialog_writes_file_and_cancel_is_reported(qapp, window, tmp_path, monkeypatch):
    wait_for(qapp, lambda: formula(qapp, window) == "H2O")
    out = tmp_path / "saved-water.xyz"
    monkeypatch.setattr(QFileDialog, "getSaveFileName", staticmethod(lambda *a, **k: (str(out), "")))
    js(qapp, window, "document.querySelector('#save-format').value = 'xyz'; document.querySelector('#btn-save').click()")
    wait_for(qapp, lambda: out.exists() and out.stat().st_size > 0)
    assert out.read_text().splitlines()[0] == "3"
    wait_for(qapp, lambda: "Saved" in (js(qapp, window, "document.querySelector('#message').textContent") or ""))
    wait_for(qapp, lambda: "saved-water" in window.windowTitle())

    monkeypatch.setattr(QFileDialog, "getSaveFileName", staticmethod(lambda *a, **k: ("", "")))
    js(qapp, window, "document.querySelector('#btn-save').click()")
    wait_for(qapp, lambda: "Save cancelled" in (js(qapp, window, "document.querySelector('#message').textContent") or ""))


def test_close_asks_about_unsaved_changes(qapp, window, monkeypatch):
    wait_for(qapp, lambda: formula(qapp, window) == "H2O")
    js(qapp, window, "document.querySelector('#smiles').value = 'C'; document.querySelector('#smiles-form').requestSubmit()")
    wait_for(qapp, lambda: formula(qapp, window) == "CH6O")
    asked = []
    monkeypatch.setattr(QMessageBox, "question", staticmethod(lambda *a, **k: asked.append(1) or QMessageBox.StandardButton.Cancel))
    window.close()
    wait_for(qapp, lambda: asked)
    assert window.isVisible(), "Cancel keeps the window open"
    monkeypatch.setattr(QMessageBox, "question", staticmethod(lambda *a, **k: QMessageBox.StandardButton.Discard))
    window.close()
    wait_for(qapp, lambda: not window.isVisible())


def _insert(qapp, win, smiles, formula_after):
    js(qapp, win, f"document.querySelector('#smiles').value = {smiles!r}; document.querySelector('#smiles-form').requestSubmit()")
    wait_for(qapp, lambda: formula(qapp, win) == formula_after)


def _move_to(qapp, win, client_xy):
    from PySide6.QtCore import QPoint

    target = win.view.focusProxy() or win.view
    QTest.mouseMove(target, QPoint(round(client_xy["x"]), round(client_xy["y"])))
    for _ in range(20):
        qapp.processEvents()
        time.sleep(0.01)


def test_hotkeys_on_atoms_and_bonds_in_the_window(qapp, window):
    wait_for(qapp, lambda: formula(qapp, window) == "H2O")
    js(qapp, window, "document.querySelector('#btn-clear').click()")
    _insert(qapp, window, "CC", "C2H6")
    target = window.view.focusProxy() or window.view
    # Hover the second carbon and press 3: it becomes a phenyl ring (toluene).
    xy = json.loads(js(qapp, window, "(() => { const w = window.__minimol; return JSON.stringify(w.viewer.modelToClient(w.M.pos(w.model().atoms[1]))); })()"))
    _move_to(qapp, window, xy)
    wait_for(qapp, lambda: js(qapp, window, "JSON.stringify(window.__minimol.state.target)") == '{"kind":"atom","index":1}')
    QTest.keyClick(target, Qt.Key.Key_3)
    wait_for(qapp, lambda: formula(qapp, window) == "C7H8")
    # Hover the CH3-C bond... find a C-C single bond whose midpoint picks it, press 2 is refused only for C-H;
    # use the ring: hover an aromatic C-C bond and press a (fuse benzene).
    xy = json.loads(js(qapp, window, """(() => {
        const w = window.__minimol, m = w.model();
        const deg = (i) => m.bonds.filter((b) => b.a === i || b.b === i).length;
        const hs = (i) => m.bonds.filter((b) => (b.a === i && m.atoms[b.b].el === 'H') || (b.b === i && m.atoms[b.a].el === 'H')).length;
        for (let k = 0; k < m.bonds.length; k++) {
          const b = m.bonds[k];
          if (m.atoms[b.a].el !== 'C' || m.atoms[b.b].el !== 'C' || hs(b.a) !== 1 || hs(b.b) !== 1 || deg(b.a) !== 3 || deg(b.b) !== 3) continue;
          const p = w.viewer.modelToClient(w.M.pos(m.atoms[b.a])), q = w.viewer.modelToClient(w.M.pos(m.atoms[b.b]));
          const mid = { x: (p.x + q.x) / 2, y: (p.y + q.y) / 2 };
          const pick = w.viewer.pick(mid.x, mid.y);
          if (pick && pick.kind === 'bond' && pick.index === k) return JSON.stringify(mid);
        }
        return "null";
    })()"""))
    assert xy is not None
    _move_to(qapp, window, xy)
    wait_for(qapp, lambda: '"kind":"bond"' in (js(qapp, window, "JSON.stringify(window.__minimol.state.target)") or ""))
    QTest.keyClick(target, Qt.Key.Key_A)
    wait_for(qapp, lambda: formula(qapp, window) == "C11H10")


def test_copy_puts_chemical_formats_on_the_system_clipboard(qapp, window):
    wait_for(qapp, lambda: formula(qapp, window) == "H2O")
    js(qapp, window, "document.querySelector('#btn-clear').click()")
    _insert(qapp, window, "CCO", "C2H6O")
    window.view.setFocus()
    QTest.keyClick(window.view.focusProxy() or window.view, Qt.Key.Key_C, Qt.KeyboardModifier.ControlModifier)
    wait_for(qapp, lambda: QApplication.clipboard().text() == "CCO")
    formats = QApplication.clipboard().mimeData().formats()
    for fmt in ("chemical/x-mdl-molfile", "chemical/x-cdxml", clipboard.MDLCT_WINDOWS_FORMAT):
        assert fmt in formats


def test_paste_reads_chemdraw_style_mdlct_from_the_system_clipboard(qapp, window):
    from PySide6.QtCore import QByteArray, QMimeData
    from rdkit import Chem

    wait_for(qapp, lambda: formula(qapp, window) == "H2O")
    data = QMimeData()
    data.setData(clipboard.MDLCT_WINDOWS_FORMAT, QByteArray(clipboard.mdlct_encode(Chem.MolToMolBlock(Chem.MolFromSmiles("CC#N")))))
    data.setData("image/png", QByteArray(b"\x89PNG"))
    QApplication.clipboard().setMimeData(data)
    window.view.setFocus()
    QTest.keyClick(window.view.focusProxy() or window.view, Qt.Key.Key_V, Qt.KeyboardModifier.ControlModifier)
    wait_for(qapp, lambda: formula(qapp, window) == "C2H5NO")  # water + acetonitrile
    assert "MDL CT" in (js(qapp, window, "document.querySelector('#message').textContent") or "")


def test_unreadable_clipboard_reports_formats(qapp, window):
    from PySide6.QtCore import QByteArray, QMimeData

    wait_for(qapp, lambda: formula(qapp, window) == "H2O")
    data = QMimeData()
    data.setData('application/x-qt-windows-mime;value="ChemDraw Interchange Format"', QByteArray(clipboard.CDX_MAGIC + b"\x00" * 64))
    QApplication.clipboard().setMimeData(data)
    window.view.setFocus()
    QTest.keyClick(window.view.focusProxy() or window.view, Qt.Key.Key_V, Qt.KeyboardModifier.ControlModifier)
    expected = "binary CDX" if not clipboard.cdx_supported() else "Nothing on the clipboard"
    wait_for(qapp, lambda: expected in (js(qapp, window, "document.querySelector('#message').textContent") or ""))


# ------------------------------------------- restart with software rendering


def close_window(win):
    win._close_confirmed = True
    win.close()
    win.deleteLater()
    QApplication.sendPostedEvents(None, QEvent.Type.DeferredDelete)


def page(script):
    """A stand-in for the editor page whose script sets the startup state."""
    import urllib.parse

    return "data:text/html," + urllib.parse.quote(f"<!doctype html><title>t</title><script>{script}</script>")


def test_window_without_webgl_closes_to_restart_with_software_rendering(qapp):
    win = desktop.create_window(page("window.minimolWebGLFailed = true;"), retry_with_software_gl=True)
    win.show()
    try:
        wait_for(qapp, lambda: win.wants_software_gl, timeout=30)
        assert win.page is None  # the window closed itself
    finally:
        close_window(win)


def test_window_already_in_software_mode_does_not_restart_again(qapp, tmp_path, monkeypatch):
    report = tmp_path / "ready.json"
    monkeypatch.setenv(desktop.READY_REPORT_ENV, str(report))
    win = desktop.create_window(page("window.minimolWebGLFailed = true;"), software_gl=True)
    win.show()
    try:
        wait_for(qapp, report.exists, timeout=30)
        result = json.loads(report.read_text())
        assert (result["state"], result["software_gl"], result["pid"]) == ("no-webgl", True, os.getpid())
        assert not win.wants_software_gl
    finally:
        close_window(win)


def test_ready_report_once_the_real_editor_has_started(qapp, tmp_path, monkeypatch):
    report = tmp_path / "ready.json"
    monkeypatch.setenv(desktop.READY_REPORT_ENV, str(report))
    server = make_server("127.0.0.1", 0, None)
    threading.Thread(target=server.serve_forever, daemon=True).start()
    win = desktop.create_window(server.url, software_gl=True, retry_with_software_gl=False)
    win.show()
    try:
        wait_for(qapp, report.exists, timeout=60)
        assert json.loads(report.read_text())["state"] == "ready"
    finally:
        close_window(win)
        server.shutdown()
        server.server_close()


def test_program_started_normally_opens_the_editor(tmp_path):
    """The whole program in its own process; without a usable GPU it restarts itself once."""
    import subprocess
    import sys

    molecule = tmp_path / "water.xyz"
    molecule.write_text("3\nwater\nO 0 0 0\nH 0.96 0 0\nH -0.24 0.93 0\n")
    report = tmp_path / "ready.json"
    env = {k: v for k, v in os.environ.items() if k not in ("MINIMOL_SOFTWARE_GL", "QTWEBENGINE_CHROMIUM_FLAGS")}
    env.update({desktop.READY_REPORT_ENV: str(report), "QT_QPA_PLATFORM": os.environ.get("QT_QPA_PLATFORM", "offscreen")})
    done = subprocess.run([sys.executable, "-m", "minimol", str(molecule)], env=env, timeout=180, capture_output=True, text=True)
    assert done.returncode == 0, done.stderr[-2000:]
    result = json.loads(report.read_text())
    assert result["state"] == "ready", result
    assert result["args"][-1] == str(molecule)  # the restarted copy (if any) got the file too
    print("restarted with software rendering:", result["software_gl"])

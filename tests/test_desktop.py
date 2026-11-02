"""Smoke tests for the desktop window (Qt WebEngine), run with Qt's offscreen platform."""

import os
import threading
import time

import pytest

pytest.importorskip("PySide6.QtWebEngineWidgets", reason="PySide6 with Qt WebEngine not available", exc_type=ImportError)
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PySide6.QtCore import QEvent, Qt  # noqa: E402
from PySide6.QtTest import QTest  # noqa: E402
from PySide6.QtWidgets import QApplication, QFileDialog, QMessageBox  # noqa: E402

from minimol import desktop  # noqa: E402
from minimol.server import make_server  # noqa: E402


@pytest.fixture(scope="module")
def qapp():
    # The test machine may have no GPU; allow WebGL on a software renderer.
    desktop._prepare_environment(software_gl=True)
    return QApplication.instance() or QApplication([])


@pytest.fixture()
def window(qapp, tmp_path):
    water = tmp_path / "water.xyz"
    water.write_text("3\nwater\nO 0 0 0\nH 0.96 0 0\nH -0.24 0.93 0\n")
    server = make_server("127.0.0.1", 0, str(water))
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

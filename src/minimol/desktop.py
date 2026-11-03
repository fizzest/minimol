"""minimol as a desktop application: a native window (Qt WebEngine) around the
same editor the browser version shows.

The local server runs on a background thread inside this process and is only
reachable from this machine. Keyboard shortcuts (Ctrl+R optimize, Ctrl+Z
undo, ...) are handled by the editor itself, so they behave the same in the
window and in a browser; inside this window nothing is bound to "reload".
Saving goes through a native save dialog, and copy/paste use the system
clipboard (with all the formats ChemDraw may provide; see clipboard.py).

Without a usable GPU (virtual machines, some remote desktops, blocklisted
drivers) Chromium refuses WebGL and the 3D view cannot start. The window then
restarts minimol once with software rendering (the same as --software-gl).
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
import threading

from .chem import ChemError
from .server import STATIC_DIR, make_server

ICON = STATIC_DIR / "icon.png"
# For tests of the installed program: once the editor has started (or failed
# to), write {"state", "software_gl", "pid", "args"} as JSON to this path and quit.
READY_REPORT_ENV = "MINIMOL_READY_REPORT"
# Whether the editor in the page has started: 'ready', 'no-webgl' or 'loading'.
PROBE_JS = "window.__minimol ? 'ready' : (window.minimolWebGLFailed ? 'no-webgl' : 'loading')"


def available() -> tuple[bool, str]:
    """Whether the Qt WebEngine window can be used, and why not."""
    try:
        from PySide6 import QtWebEngineWidgets  # noqa: F401
    except ImportError as exc:
        return False, str(exc)
    return True, ""


def software_gl_requested(flag: bool = False) -> bool:
    return flag or os.environ.get("MINIMOL_SOFTWARE_GL") == "1"


def software_gl_command(initial_file: str | None, port: int) -> list[str]:
    """The command that starts this program again with software rendering."""
    if getattr(sys, "frozen", False):  # the packaged program (PyInstaller)
        cmd = [sys.executable]
    else:
        cmd = [sys.executable, "-m", "minimol"]
    cmd.append("--software-gl")
    if port:
        cmd += ["--port", str(port)]
    if initial_file:
        cmd += ["--", initial_file]
    return cmd


def _start_again(cmd: list[str]) -> int:
    if sys.platform == "win32":
        # A windowed program: start the new copy and let this one exit.
        subprocess.Popen(cmd, close_fds=True)
        return 0
    # Elsewhere minimol may have been started from a terminal; stay attached to it.
    try:
        return subprocess.call(cmd)
    except KeyboardInterrupt:
        return 130


def _write_ready_report(path: str, state: str, software_gl: bool) -> None:
    tmp = f"{path}.tmp"
    with open(tmp, "w", encoding="utf-8") as fh:
        json.dump({"state": state, "software_gl": software_gl, "pid": os.getpid(), "args": sys.argv[1:]}, fh)
    os.replace(tmp, path)  # readers never see a half-written file


def _prepare_environment(software_gl: bool = False) -> None:
    """Set Qt WebEngine options; must run before the QApplication is created."""
    if sys.platform == "win32":
        # Own taskbar entry and icon (instead of being grouped as "Python").
        try:
            import ctypes

            ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID("minimol.minimol")
        except (AttributeError, OSError):
            pass
    # Chromium refuses to start its sandbox as root on Linux; without this
    # the window would stay blank. Normal user accounts keep the sandbox.
    if sys.platform.startswith("linux") and hasattr(os, "geteuid") and os.geteuid() == 0:
        os.environ.setdefault("QTWEBENGINE_DISABLE_SANDBOX", "1")
    # Without a usable GPU, Chromium blocks WebGL and the 3D view cannot
    # start. Ignoring the GPU blocklist lets it use a software renderer
    # (verified with Mesa llvmpipe on Linux).
    if software_gl_requested(software_gl):
        flags = os.environ.get("QTWEBENGINE_CHROMIUM_FLAGS", "").split()
        wanted = ["--ignore-gpu-blocklist"]
        if sys.platform.startswith("linux") and not os.path.exists("/dev/dri"):
            # No GPU device (e.g. WSL): Chromium's GPU compositor then tries to
            # share frames through dma_buf and the window stays blank.
            wanted.append("--disable-gpu-compositing")
        os.environ["QTWEBENGINE_CHROMIUM_FLAGS"] = " ".join(flags + [f for f in wanted if f not in flags])


def make_clipboard_bridge():
    """An object giving the server thread access to the system clipboard.

    Qt only allows clipboard access from the GUI thread, while the server
    answers requests on its own thread. The bridge sends each request to the
    GUI thread through a queued signal and waits for the answer. It must be
    created on the GUI thread, after the QApplication.
    """
    from PySide6.QtCore import QByteArray, QMimeData, QObject, Signal
    from PySide6.QtGui import QGuiApplication

    class ClipboardBridge(QObject):
        _request = Signal(object)

        def __init__(self):
            super().__init__()
            self._request.connect(self._handle)

        def _handle(self, job):
            try:
                clip = QGuiApplication.clipboard()
                if job["kind"] == "set":
                    data = QMimeData()
                    for fmt, raw in job["data"].items():
                        if fmt == "text/plain":
                            data.setText(raw.decode("utf-8"))
                        else:
                            data.setData(fmt, QByteArray(raw))
                    clip.setMimeData(data)
                    job["result"] = None
                else:
                    current = clip.mimeData()
                    job["result"] = {f: current.data(f).data() for f in current.formats()} if current else {}
            except Exception as exc:  # reported to the caller
                job["error"] = exc
            finally:
                job["done"].set()

        def _call(self, kind, data=None):
            job = {"kind": kind, "data": data, "done": threading.Event()}
            self._request.emit(job)
            if not job["done"].wait(10):
                raise ChemError("The system clipboard did not respond")
            if "error" in job:
                raise ChemError(f"Clipboard error: {job['error']}")
            return job["result"]

        def set(self, formats: dict[str, bytes]) -> None:
            self._call("set", formats)

        def get(self) -> dict[str, bytes]:
            return self._call("get")

    return ClipboardBridge()


def create_window(url: str, software_gl: bool = False, retry_with_software_gl: bool = False):
    """Build the main window for the editor at ``url`` (QApplication must exist).

    ``software_gl`` says whether this process renders in software;
    with ``retry_with_software_gl`` the window closes itself, setting
    ``wants_software_gl``, when the 3D view cannot start.
    """
    from PySide6.QtCore import Qt, QTimer, QUrl
    from PySide6.QtWebEngineCore import QWebEngineDownloadRequest, QWebEnginePage, QWebEngineProfile
    from PySide6.QtWebEngineWidgets import QWebEngineView
    from PySide6.QtWidgets import QFileDialog, QMainWindow, QMessageBox

    class MainWindow(QMainWindow):
        def __init__(self):
            super().__init__()
            self.setWindowTitle("minimol")
            self.resize(1360, 860)
            # A private in-memory profile per window: a shared on-disk profile
            # breaks when two minimol windows run at the same time.
            self.profile = QWebEngineProfile(self)
            self.page = QWebEnginePage(self.profile, self)
            self.view = QWebEngineView(self)
            self.view.setPage(self.page)
            # Right-click belongs to the editor (delete), not to a browser menu.
            self.view.setContextMenuPolicy(Qt.ContextMenuPolicy.NoContextMenu)
            self.setCentralWidget(self.view)
            self.view.titleChanged.connect(self.setWindowTitle)
            self.profile.downloadRequested.connect(self._save)
            self.page.loadFinished.connect(self._loaded)
            self._close_confirmed = False
            self._last_dir = os.path.expanduser("~")
            self._software_gl = software_gl
            self._retry = retry_with_software_gl
            self._report = os.environ.get(READY_REPORT_ENV)
            self._probing = False
            self.wants_software_gl = False
            self.view.load(QUrl(url))
            self.view.setFocus()

        def js(self, code: str, callback=None):
            if self.page is None:
                return
            if callback is None:
                self.page.runJavaScript(code, 0)
            else:
                self.page.runJavaScript(code, 0, callback)

        def _loaded(self, ok: bool):
            if ok:
                self.js("window.minimolDesktop = true;")
                if (self._retry or self._report) and not self._probing:
                    self._probing = True
                    self._probe(tries_left=120)
            self.view.setFocus()

        # ------------------------------------------- did the editor start?

        def _probe(self, tries_left: int):
            self.js(PROBE_JS, lambda state: self._probed(state, tries_left))

        def _probed(self, state, tries_left: int):
            if self.page is None:
                return
            if state not in ("ready", "no-webgl"):
                if tries_left > 0:
                    QTimer.singleShot(250, lambda: self._probe(tries_left - 1))
                    return
                state = "timed out"
            if state == "no-webgl" and self._retry:
                self.wants_software_gl = True  # run() starts minimol again
                self._close_confirmed = True
                self.close()
                return
            if self._report:
                _write_ready_report(self._report, state, self._software_gl)
                self._close_confirmed = True
                self.close()

        # ------------------------------------------------------------ saving

        def _save(self, download: QWebEngineDownloadRequest):
            start = os.path.join(self._last_dir, download.downloadFileName())
            path, _ = QFileDialog.getSaveFileName(self, "Save molecule", start)
            if not path:
                download.cancel()
                self.js("window.minimolSaveFailed && window.minimolSaveFailed('')")
                return
            self._last_dir = os.path.dirname(path)
            download.setDownloadDirectory(os.path.dirname(path))
            download.setDownloadFileName(os.path.basename(path))

            def finished():
                if not download.isFinished():
                    return
                if download.state() == QWebEngineDownloadRequest.DownloadState.DownloadCompleted:
                    self.js(f"window.minimolSaved && window.minimolSaved({json.dumps(path)})")
                else:
                    reason = download.interruptReasonString() or "the file could not be written"
                    self.js(f"window.minimolSaveFailed && window.minimolSaveFailed({json.dumps(reason)})")

            download.isFinishedChanged.connect(finished)
            download.accept()

        # ----------------------------------------------------------- closing

        def release_web(self):
            """Delete the page before its profile (Qt WebEngine requires this order)."""
            if self.page is None:
                return
            self.view.setPage(QWebEnginePage(self.view))  # placeholder on the default profile
            self.page.deleteLater()
            self.profile.deleteLater()
            self.page = None

        def closeEvent(self, event):
            if self._close_confirmed:
                self.release_web()
                event.accept()
                return
            event.ignore()
            self.js(
                "window.minimolHasUnsavedChanges ? window.minimolHasUnsavedChanges() : false",
                self._close_answer,
            )

        def _close_answer(self, unsaved):
            if unsaved:
                answer = QMessageBox.question(
                    self,
                    "minimol",
                    "The molecule has unsaved changes. Close anyway?",
                    QMessageBox.StandardButton.Discard | QMessageBox.StandardButton.Cancel,
                    QMessageBox.StandardButton.Cancel,
                )
                if answer != QMessageBox.StandardButton.Discard:
                    return
            self._close_confirmed = True
            self.close()

    return MainWindow()


def run(initial_file: str | None = None, port: int = 0, software_gl: bool = False) -> int:
    """Start the server and the desktop window; returns the exit code."""
    software_gl = software_gl_requested(software_gl)
    _prepare_environment(software_gl)
    from PySide6.QtCore import QEvent
    from PySide6.QtWidgets import QApplication

    server = make_server("127.0.0.1", port, initial_file)
    thread = threading.Thread(target=server.serve_forever, name="minimol-server", daemon=True)
    thread.start()
    from PySide6.QtGui import QIcon

    app = QApplication.instance() or QApplication(sys.argv)
    app.setApplicationName("minimol")
    app.setApplicationDisplayName("minimol")
    if ICON.exists():
        app.setWindowIcon(QIcon(str(ICON)))
    server.clipboard = make_clipboard_bridge()
    window = create_window(server.url, software_gl=software_gl, retry_with_software_gl=not software_gl)
    window.show()
    try:
        code = app.exec()
    finally:
        window.release_web()
        QApplication.sendPostedEvents(None, QEvent.Type.DeferredDelete)
        server.shutdown()
        server.server_close()
    if window.wants_software_gl:
        # Only once: the new copy runs with --software-gl and never retries.
        return _start_again(software_gl_command(initial_file, port))
    return code

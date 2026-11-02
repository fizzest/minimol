"""minimol as a desktop application: a native window (Qt WebEngine) around the
same editor the browser version shows.

The local server runs on a background thread inside this process and is only
reachable from this machine. Keyboard shortcuts (Ctrl+R optimize, Ctrl+Z
undo, ...) are handled by the editor itself, so they behave the same in the
window and in a browser; inside this window nothing is bound to "reload".
Saving goes through a native save dialog.
"""

from __future__ import annotations

import json
import os
import sys
import threading

from .server import make_server


def available() -> tuple[bool, str]:
    """Whether the Qt WebEngine window can be used, and why not."""
    try:
        from PySide6 import QtWebEngineWidgets  # noqa: F401
    except ImportError as exc:
        return False, str(exc)
    return True, ""


def _prepare_environment(software_gl: bool = False) -> None:
    """Set Qt WebEngine options; must run before the QApplication is created."""
    # Chromium refuses to start its sandbox as root on Linux; without this
    # the window would stay blank. Normal user accounts keep the sandbox.
    if sys.platform.startswith("linux") and hasattr(os, "geteuid") and os.geteuid() == 0:
        os.environ.setdefault("QTWEBENGINE_DISABLE_SANDBOX", "1")
    # Without a usable GPU, Chromium blocks WebGL and the 3D view cannot
    # start. Ignoring the GPU blocklist lets it use a software renderer
    # (verified with Mesa llvmpipe on Linux).
    if software_gl or os.environ.get("MINIMOL_SOFTWARE_GL") == "1":
        flags = os.environ.get("QTWEBENGINE_CHROMIUM_FLAGS", "")
        if "--ignore-gpu-blocklist" not in flags:
            os.environ["QTWEBENGINE_CHROMIUM_FLAGS"] = (flags + " --ignore-gpu-blocklist").strip()


def create_window(url: str):
    """Build the main window for the editor at ``url`` (QApplication must exist)."""
    from PySide6.QtCore import Qt, QUrl
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
            self.view.setFocus()

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
    _prepare_environment(software_gl)
    from PySide6.QtCore import QEvent
    from PySide6.QtWidgets import QApplication

    server = make_server("127.0.0.1", port, initial_file)
    thread = threading.Thread(target=server.serve_forever, name="minimol-server", daemon=True)
    thread.start()
    app = QApplication.instance() or QApplication(sys.argv)
    app.setApplicationName("minimol")
    app.setApplicationDisplayName("minimol")
    window = create_window(server.url)
    window.show()
    try:
        return app.exec()
    finally:
        window.release_web()
        QApplication.sendPostedEvents(None, QEvent.Type.DeferredDelete)
        server.shutdown()
        server.server_close()

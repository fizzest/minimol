"""The desktop app's restart with software rendering (no Qt needed for these)."""

import sys

from minimol import __main__ as cli
from minimol import desktop


def test_installed_program_restarts_itself(monkeypatch):
    exe = r"C:\Users\me\AppData\Local\Programs\minimol\minimol.exe"
    monkeypatch.setattr(desktop.sys, "frozen", True, raising=False)
    monkeypatch.setattr(desktop.sys, "executable", exe)
    assert desktop.software_gl_command(r"C:\data\my molecule.xyz", 0) == [exe, "--software-gl", "--", r"C:\data\my molecule.xyz"]
    assert desktop.software_gl_command(None, 0) == [exe, "--software-gl"]


def test_restart_command_keeps_the_file_and_port(monkeypatch, tmp_path):
    monkeypatch.delattr(desktop.sys, "frozen", raising=False)
    molecule = tmp_path / "water.xyz"
    molecule.write_text("3\nwater\nO 0 0 0\nH 0.96 0 0\nH -0.24 0.93 0\n")
    cmd = desktop.software_gl_command(str(molecule), 8123)
    assert cmd[:3] == [sys.executable, "-m", "minimol"]
    # The new copy's command line means the same thing to minimol's parser.
    seen = {}
    monkeypatch.setattr(desktop, "available", lambda: (True, ""))
    monkeypatch.setattr(desktop, "run", lambda initial, port, software_gl: seen.update(file=initial, port=port, software_gl=software_gl) or 0)
    assert cli.main(cmd[3:]) == 0
    assert seen == {"file": str(molecule), "port": 8123, "software_gl": True}


def test_restart_on_windows_detaches_and_elsewhere_waits(monkeypatch):
    calls = []
    monkeypatch.setattr(desktop.subprocess, "call", lambda cmd: calls.append(("wait", cmd)) or 3)
    monkeypatch.setattr(desktop.subprocess, "Popen", lambda cmd, **kw: calls.append(("detach", cmd)))
    monkeypatch.setattr(desktop.sys, "platform", "linux")
    assert desktop._start_again(["minimol"]) == 3  # the new copy's exit code
    monkeypatch.setattr(desktop.sys, "platform", "win32")
    assert desktop._start_again(["minimol.exe"]) == 0
    assert calls == [("wait", ["minimol"]), ("detach", ["minimol.exe"])]


def test_software_gl_can_also_come_from_the_environment(monkeypatch):
    monkeypatch.delenv("MINIMOL_SOFTWARE_GL", raising=False)
    assert not desktop.software_gl_requested(False)
    assert desktop.software_gl_requested(True)
    monkeypatch.setenv("MINIMOL_SOFTWARE_GL", "1")
    assert desktop.software_gl_requested(False)

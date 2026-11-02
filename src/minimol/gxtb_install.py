"""Download g-xTB for this computer, on request ("Install g-xTB..." in the app).

g-xTB is GPL-3.0 software published by the Grimme group; minimol does not
bundle it. When the user asks, this module downloads the official archive for
the operating system from the g-xTB repository, pinned to one commit, checks
it against the SHA-256 published next to it in that repository (copied here),
and unpacks it into minimol's own data folder:

* Windows: %LOCALAPPDATA%\\minimol\\g-xtb
* macOS: ~/Library/Application Support/minimol/g-xtb
* Linux: $XDG_DATA_HOME/minimol/g-xtb (default ~/.local/share/minimol/g-xtb)

gxtb.find_binary() then finds it there.
"""

from __future__ import annotations

import hashlib
import os
import platform
import shutil
import sys
import tarfile
import tempfile
import urllib.request
import zipfile
from pathlib import Path

from .chem import ChemError

REPOSITORY = "https://github.com/grimme-lab/g-xtb"
COMMIT = "23f81c2db64fa4272c25b86f8776ce862aeb59e2"  # "Version 2.0.1 - small bugfixes", 2026-05-14
BASE_URL = f"https://raw.githubusercontent.com/grimme-lab/g-xtb/{COMMIT}/binaries/"
# (sys.platform, machine) -> (archive name, SHA-256 from the repository's .sha256 file)
ARCHIVES = {
    ("win32", "amd64"): ("xtb-6.7.1-gxtb-140526-windows-x86_64.zip", "fa6d6491b38d895196e2312c6ce34b74e60cd974f9c01da6c3ee567b3ca41830"),
    ("linux", "x86_64"): ("xtb-6.7.1-gxtb-140526-linux-x86_64.tar.xz", "cdfb81164f8d3300fc451dbbf2414528b533e71f7a30ec824b564cd30d6d9e5c"),
    ("darwin", "arm64"): ("xtb-6.7.1-gxtb-140526-macos-arm64.tar.gz", "035a37f241a24067eb854177ea0dd3e1d1b48aca2155cf934aae315351e098af"),
    ("darwin", "x86_64"): ("xtb-6.7.1-gxtb-140526-macos-x86_64.tar.gz", "3f5f7ee80a66bce911cd51dcca894fa087e7039b51e8c2e036ba74208a8e3b62"),
}
MAX_DOWNLOAD = 200 * 1024 * 1024  # the archives are 35-45 MB


def data_dir() -> Path:
    """minimol's per-user data folder."""
    if sys.platform == "win32":
        root = os.environ.get("LOCALAPPDATA") or os.path.expanduser(r"~\AppData\Local")
        return Path(root) / "minimol"
    if sys.platform == "darwin":
        return Path.home() / "Library" / "Application Support" / "minimol"
    root = os.environ.get("XDG_DATA_HOME") or os.path.expanduser("~/.local/share")
    return Path(root) / "minimol"


def install_dir() -> Path:
    return data_dir() / "g-xtb"


def _binary_name() -> str:
    return "xtb.exe" if sys.platform == "win32" else "xtb"


def installed_binary() -> Path | None:
    """The xtb binary minimol installed earlier, if any."""
    root = install_dir()
    if not root.is_dir():
        return None
    for path in sorted(root.rglob(_binary_name())):
        if path.parent.name == "bin" and path.is_file():
            return path
    return None


def archive_for_this_computer() -> tuple[str, str] | None:
    machine = platform.machine().lower()
    machine = {"x86_64": "amd64" if sys.platform == "win32" else "x86_64", "aarch64": "arm64"}.get(machine, machine)
    return ARCHIVES.get((sys.platform, machine))


def _download(url: str, dest: Path, expected_sha256: str, opener=urllib.request.urlopen) -> None:
    digest = hashlib.sha256()
    size = 0
    try:
        with opener(url, timeout=60) as response, open(dest, "wb") as out:
            while True:
                chunk = response.read(1 << 16)
                if not chunk:
                    break
                size += len(chunk)
                if size > MAX_DOWNLOAD:
                    raise ChemError("The g-xTB download is unexpectedly large; stopped.")
                digest.update(chunk)
                out.write(chunk)
    except ChemError:
        raise
    except OSError as exc:
        raise ChemError(f"Could not download g-xTB: {exc}") from None
    if digest.hexdigest() != expected_sha256:
        raise ChemError("The downloaded g-xTB archive does not match its published checksum; nothing was installed.")


def _safe_target(root: Path, name: str) -> Path:
    target = (root / name).resolve()
    if target != root.resolve() and root.resolve() not in target.parents:
        raise ChemError(f"The g-xTB archive contains an unsafe path: {name}")
    return target


def _extract(archive: Path, dest: Path) -> None:
    if archive.name.endswith(".zip"):
        with zipfile.ZipFile(archive) as z:
            for info in z.infolist():
                target = _safe_target(dest, info.filename)
                if info.is_dir():
                    target.mkdir(parents=True, exist_ok=True)
                    continue
                target.parent.mkdir(parents=True, exist_ok=True)
                with z.open(info) as src, open(target, "wb") as out:
                    shutil.copyfileobj(src, out)
        return
    with tarfile.open(archive) as t:
        for member in t.getmembers():
            target = _safe_target(dest, member.name)
            if member.isdir():
                target.mkdir(parents=True, exist_ok=True)
            elif member.isfile():
                target.parent.mkdir(parents=True, exist_ok=True)
                with t.extractfile(member) as src, open(target, "wb") as out:
                    shutil.copyfileobj(src, out)
                if member.mode & 0o111:
                    target.chmod(0o755)
            # Links and special files are skipped: the archive does not need them.


def install(opener=urllib.request.urlopen) -> Path:
    """Download, verify and unpack g-xTB; returns the path of its xtb binary."""
    choice = archive_for_this_computer()
    if choice is None:
        raise ChemError(f"There is no g-xTB build for this computer ({sys.platform}, {platform.machine()}).")
    name, sha256 = choice
    target = install_dir()
    target.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="minimol-gxtb-", dir=target.parent) as work:
        work = Path(work)
        archive = work / name
        _download(BASE_URL + name, archive, sha256, opener)
        unpacked = work / "unpacked"
        unpacked.mkdir()
        _extract(archive, unpacked)
        found = [p for p in unpacked.rglob(_binary_name()) if p.parent.name == "bin"]
        if not found:
            raise ChemError("The g-xTB archive does not contain the expected bin/xtb program.")
        if target.exists():
            shutil.rmtree(target)
        shutil.move(str(unpacked), str(target))
    binary = installed_binary()
    if binary is None:  # pragma: no cover - would mean the move failed silently
        raise ChemError("g-xTB was unpacked but its program could not be found.")
    return binary

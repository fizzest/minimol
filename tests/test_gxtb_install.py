import hashlib
import io
import os
import stat
import tarfile
import zipfile

import pytest

from minimol import gxtb, gxtb_install
from minimol.chem import ChemError


class FakeResponse(io.BytesIO):
    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()


def fake_opener(payload: bytes, seen: list):
    def opener(url, timeout=None):
        seen.append(url)
        return FakeResponse(payload)

    return opener


def tar_xz(files: dict[str, bytes], mode=0o755) -> bytes:
    buf = io.BytesIO()
    with tarfile.open(fileobj=buf, mode="w:xz") as t:
        for name, data in files.items():
            info = tarfile.TarInfo(name)
            info.size = len(data)
            info.mode = mode
            t.addfile(info, io.BytesIO(data))
    return buf.getvalue()


@pytest.fixture()
def linux_target(tmp_path, monkeypatch):
    """Pretend to be Linux x86-64 with a private data folder, and use a fake archive."""
    monkeypatch.setattr(gxtb_install.sys, "platform", "linux")
    monkeypatch.setattr(gxtb_install.platform, "machine", lambda: "x86_64")
    monkeypatch.setenv("XDG_DATA_HOME", str(tmp_path / "data"))
    return tmp_path


def use_archive(monkeypatch, payload: bytes):
    name = gxtb_install.ARCHIVES[("linux", "x86_64")][0]
    monkeypatch.setitem(gxtb_install.ARCHIVES, ("linux", "x86_64"), (name, hashlib.sha256(payload).hexdigest()))


def test_pinned_archives_match_the_published_checksums():
    # The four archives and checksums of g-xTB 2.0.1 (repository .sha256 files at the pinned commit).
    assert gxtb_install.COMMIT == "23f81c2db64fa4272c25b86f8776ce862aeb59e2"
    assert gxtb_install.ARCHIVES[("win32", "amd64")] == (
        "xtb-6.7.1-gxtb-140526-windows-x86_64.zip",
        "fa6d6491b38d895196e2312c6ce34b74e60cd974f9c01da6c3ee567b3ca41830",
    )
    assert all(len(sha) == 64 for _, sha in gxtb_install.ARCHIVES.values())
    assert gxtb_install.BASE_URL.startswith("https://raw.githubusercontent.com/grimme-lab/g-xtb/" + gxtb_install.COMMIT)


def test_install_downloads_verifies_and_unpacks(linux_target, monkeypatch):
    payload = tar_xz({"xtb-6.7.1/bin/xtb": b"#!/bin/sh\nexit 0\n", "xtb-6.7.1/share/readme": b"hi"})
    use_archive(monkeypatch, payload)
    seen = []
    binary = gxtb_install.install(opener=fake_opener(payload, seen))
    assert seen == [gxtb_install.BASE_URL + "xtb-6.7.1-gxtb-140526-linux-x86_64.tar.xz"]
    assert binary == gxtb_install.installed_binary()
    assert binary.read_bytes().startswith(b"#!/bin/sh")
    if os.name != "nt":  # Windows has no executable bit
        assert binary.stat().st_mode & stat.S_IXUSR
    assert str(binary).startswith(str(linux_target / "data" / "minimol" / "g-xtb"))
    # Installing again replaces the old copy.
    assert gxtb_install.install(opener=fake_opener(payload, seen)) == binary


def test_checksum_mismatch_installs_nothing(linux_target, monkeypatch):
    good = tar_xz({"xtb-6.7.1/bin/xtb": b"x"})
    use_archive(monkeypatch, good)
    with pytest.raises(ChemError, match="checksum"):
        gxtb_install.install(opener=fake_opener(good + b"tampered", []))
    assert gxtb_install.installed_binary() is None


def test_unsafe_archive_paths_are_rejected(linux_target, monkeypatch):
    evil = tar_xz({"../../escaped": b"x", "xtb/bin/xtb": b"x"})
    use_archive(monkeypatch, evil)
    with pytest.raises(ChemError, match="unsafe path"):
        gxtb_install.install(opener=fake_opener(evil, []))
    assert not (linux_target / "escaped").exists()


def test_zip_archives_unpack(tmp_path):
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w") as z:
        z.writestr("xtb-bleed-windows/bin/xtb.exe", b"MZ")
        z.writestr("xtb-bleed-windows/bin/libgfortran-5.dll", b"MZ")
    archive = tmp_path / "a.zip"
    archive.write_bytes(buf.getvalue())
    gxtb_install._extract(archive, tmp_path / "out")
    assert (tmp_path / "out" / "xtb-bleed-windows" / "bin" / "xtb.exe").read_bytes() == b"MZ"


def test_unsupported_computer_is_reported(monkeypatch):
    monkeypatch.setattr(gxtb_install.sys, "platform", "sunos5")
    with pytest.raises(ChemError, match="no g-xTB build"):
        gxtb_install.install(opener=fake_opener(b"", []))


@pytest.mark.skipif(os.name == "nt", reason="uses a shell script as a stand-in for g-xTB")
def test_find_binary_uses_the_installed_copy(linux_target, monkeypatch):
    script = b"#!/bin/sh\nif [ \"$1\" = \"--help\" ]; then echo '--gxtb  use g-xTB'; fi\n"
    payload = tar_xz({"pkg/bin/xtb": script})
    use_archive(monkeypatch, payload)
    binary = gxtb_install.install(opener=fake_opener(payload, []))
    monkeypatch.setenv("PATH", str(linux_target / "empty"))
    monkeypatch.delenv(gxtb.ENV_VAR, raising=False)
    try:
        assert gxtb.find_binary(refresh=True) == str(binary)
    finally:
        gxtb._cache.clear()


def test_windows_runs_g_xtb_without_a_console_window():
    if os.name == "nt":
        assert gxtb._NO_WINDOW["creationflags"]
    else:
        assert gxtb._NO_WINDOW == {}

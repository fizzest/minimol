# PyInstaller spec for the minimol desktop program.
#
# Build (from the repository root, with minimol and pyinstaller installed):
#   pyinstaller --noconfirm --clean --distpath dist --workpath build packaging/minimol.spec
# Result: dist/minimol/ with minimol.exe (Windows) or minimol (Linux), which the
# Inno Setup script packaging/windows/minimol.iss turns into an installer.
import sys
from pathlib import Path

from PyInstaller.utils.hooks import collect_all, collect_data_files

ROOT = Path(SPECPATH).resolve().parent

datas = collect_data_files("minimol") + [(str(ROOT / "packaging" / "THIRD_PARTY_NOTICES.txt"), "."), (str(ROOT / "LICENSE"), ".")]
binaries = []
hiddenimports = ["minimol.desktop", "minimol.selftest"]
rd_datas, rd_binaries, rd_hidden = collect_all("rdkit")
datas += rd_datas
binaries += rd_binaries
hiddenimports += rd_hidden

a = Analysis(
    [str(ROOT / "packaging" / "launcher.py")],
    pathex=[str(ROOT / "src")],
    binaries=binaries,
    datas=datas,
    hiddenimports=hiddenimports,
    # Not used by minimol; keeping them out makes the program smaller.
    excludes=["tkinter", "matplotlib", "IPython", "pandas", "PIL", "pytest", "PyQt5", "PyQt6"],
)
pyz = PYZ(a.pure)
exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name="minimol",
    console=False,  # a windowed program: no console window on Windows
    icon=str(ROOT / "packaging" / "windows" / "minimol.ico") if sys.platform == "win32" else None,
)
coll = COLLECT(exe, a.binaries, a.datas, strip=False, upx=False, name="minimol")

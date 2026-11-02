"""Entry point of the packaged (PyInstaller) minimol program."""

import sys


def _log_to_file() -> None:
    # The windowed Windows program has no console: sys.stdout/stderr are None.
    # Send anything printed (e.g. a server error) to a log file instead.
    if sys.stdout is not None and sys.stderr is not None:
        return
    from minimol.gxtb_install import data_dir

    folder = data_dir()
    folder.mkdir(parents=True, exist_ok=True)
    log = open(folder / "minimol.log", "a", encoding="utf-8", buffering=1)  # noqa: SIM115 - lives as long as the program
    sys.stdout = sys.stdout or log
    sys.stderr = sys.stderr or log


if __name__ == "__main__":
    _log_to_file()
    from minimol.__main__ import main

    sys.exit(main())

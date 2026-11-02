"""Command line entry point: ``minimol [FILE] [--port N] [--no-browser]``."""

from __future__ import annotations

import argparse
import os
import sys
import webbrowser

from . import __version__
from .server import make_server


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="minimol", description="Minimal 3D molecule viewer and editor.")
    parser.add_argument("file", nargs="?", help="molecule file to open (.xyz, .mol, .sdf, .mol2, .pdb)")
    parser.add_argument("--port", type=int, default=8765, help="port to listen on (default 8765)")
    parser.add_argument("--no-browser", action="store_true", help="do not open a browser window")
    parser.add_argument("--version", action="version", version=f"minimol {__version__}")
    args = parser.parse_args(argv)

    if args.file and not os.path.isfile(args.file):
        parser.error(f"no such file: {args.file}")
    if not 0 <= args.port <= 65535:
        parser.error("port must be between 0 and 65535")

    server = make_server("127.0.0.1", args.port, os.path.abspath(args.file) if args.file else None)
    print(f"minimol {__version__} running at {server.url}  (Ctrl+C to quit)", flush=True)
    if not args.no_browser:
        webbrowser.open(server.url)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        server.server_close()
    return 0


if __name__ == "__main__":
    sys.exit(main())

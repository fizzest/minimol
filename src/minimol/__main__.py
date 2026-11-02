"""Command line entry point.

    minimol [FILE]                 open the desktop window
    minimol --browser [FILE]       use your web browser instead
    minimol --software-gl [FILE]   desktop window without a usable GPU
"""

from __future__ import annotations

import argparse
import os
import sys
import webbrowser

from . import __version__


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="minimol", description="Minimal 3D molecule viewer and editor.")
    parser.add_argument("file", nargs="?", help="molecule file to open (.xyz, .mol, .sdf, .mol2, .pdb)")
    parser.add_argument("--browser", action="store_true", help="open in your web browser instead of the desktop window")
    parser.add_argument("--port", type=int, default=None, help="port for the local server (browser mode default 8765)")
    parser.add_argument("--no-browser", action="store_true", help="browser mode: only start the server")
    parser.add_argument(
        "--software-gl",
        action="store_true",
        help="desktop window: allow software 3D rendering when the GPU is unavailable or blocked",
    )
    parser.add_argument("--version", action="version", version=f"minimol {__version__}")
    args = parser.parse_args(argv)

    if args.file and not os.path.isfile(args.file):
        parser.error(f"no such file: {args.file}")
    if args.port is not None and not 0 <= args.port <= 65535:
        parser.error("port must be between 0 and 65535")
    initial = os.path.abspath(args.file) if args.file else None

    if not args.browser:
        from . import desktop

        ok, reason = desktop.available()
        if ok:
            return desktop.run(initial, args.port or 0, args.software_gl)
        print(f"Desktop window unavailable ({reason}); opening in the browser instead.", file=sys.stderr)

    from .server import make_server

    server = make_server("127.0.0.1", 8765 if args.port is None else args.port, initial)
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

"""Command line: `pe ui` starts the app; `pe checks` reruns every model check."""

from __future__ import annotations

import argparse
import socket
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
FIRST_PORT = 8502  # 8501 is Streamlit's default, often held by another app (catalyst-engine's)


def free_port(start: int = FIRST_PORT, tries: int = 50) -> int:
    """The first port from `start` that this machine can listen on."""
    for port in range(start, start + tries):
        with socket.socket() as s:
            try:
                s.bind(("", port))
            except OSError:
                continue
            return port
    raise SystemExit(f"Ports {start}-{start + tries - 1} are all in use; pass --port.")


def main() -> None:
    parser = argparse.ArgumentParser(prog="pe", description=__doc__)
    sub = parser.add_subparsers(dest="cmd", required=True)
    ui = sub.add_parser("ui", help="start the app")
    ui.add_argument("--port", type=int, help=f"default: the first free port from {FIRST_PORT}")
    sub.add_parser("checks", help="rerun every check in experiments/")
    args = parser.parse_args()
    if args.cmd == "ui":
        app = Path(__file__).parent / "ui" / "app.py"
        port = args.port or free_port()
        sys.exit(subprocess.call([sys.executable, "-m", "streamlit", "run", str(app), "--server.port", str(port)],
                                 cwd=ROOT))
    for script in sorted((ROOT / "experiments").glob("*/check.py")):
        print(f"== {script.parent.name}")
        code = subprocess.call([sys.executable, str(script)], cwd=ROOT)
        if code:
            sys.exit(code)

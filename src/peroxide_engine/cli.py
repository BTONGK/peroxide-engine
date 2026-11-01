"""Command line: `pe ui` starts the app; `pe checks` reruns every model check."""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def main() -> None:
    parser = argparse.ArgumentParser(prog="pe", description=__doc__)
    sub = parser.add_subparsers(dest="cmd", required=True)
    ui = sub.add_parser("ui", help="start the app")
    ui.add_argument("--port", type=int, default=8501)
    sub.add_parser("checks", help="rerun every check in experiments/")
    args = parser.parse_args()
    if args.cmd == "ui":
        app = Path(__file__).parent / "ui" / "app.py"
        sys.exit(subprocess.call([sys.executable, "-m", "streamlit", "run", str(app), "--server.port", str(args.port)],
                                 cwd=ROOT))
    for script in sorted((ROOT / "experiments").glob("*/check.py")):
        print(f"== {script.parent.name}")
        code = subprocess.call([sys.executable, str(script)], cwd=ROOT)
        if code:
            sys.exit(code)

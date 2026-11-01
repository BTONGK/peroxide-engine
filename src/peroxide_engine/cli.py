"""Command line: `pe ui` starts the app; `pe checks` reruns every model check; `pe screen` prepares a computed screen
of Materials Project crystals as job folders for the UMA worker."""

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


def screen(args) -> None:
    from peroxide_engine.compute import jobs, materials

    filters = dict(elements=args.elements, exclude=args.exclude, max_elements=args.max_elements,
                   max_e_hull=args.max_hull, max_sites=args.max_sites)
    n = materials.count(**filters)
    print(f"{n:,} Materials Project crystals match")
    if args.count or n == 0:
        return
    rows = materials.search(**filters, limit=args.limit)
    out = args.out or ROOT / "data" / "screens" / f"{'-'.join(args.elements) or 'all'}-{len(rows)}"
    print(f"Cutting the faces of {len(rows):,} crystals into {out}")
    m = jobs.prepare_screen(rows, out, args.checkpoint, args.task, args.adsorbates, per_job=args.per_job)
    print(f"{m['ready']:,} ready with {m['faces']:,} faces in {m['shards']} job folders; {len(m['set_aside'])} set "
          f"aside (listed in manifest.json). Run each folder with worker.py, or sbatch run_slurm.sh on a cluster.")


def main() -> None:
    parser = argparse.ArgumentParser(prog="pe", description=__doc__)
    sub = parser.add_subparsers(dest="cmd", required=True)
    ui = sub.add_parser("ui", help="start the app")
    ui.add_argument("--port", type=int, help=f"default: the first free port from {FIRST_PORT}")
    sub.add_parser("checks", help="rerun every check in experiments/")
    sc = sub.add_parser("screen", help="select Materials Project crystals and write UMA job folders")
    sc.add_argument("--elements", nargs="*", default=["O"], help="elements every material contains (default: O)")
    sc.add_argument("--exclude", nargs="*", default=[], help="elements to leave out")
    sc.add_argument("--max-elements", type=int, help="most distinct elements per material")
    sc.add_argument("--max-hull", type=float, default=0.05, help="eV/atom above the convex hull (default 0.05)")
    sc.add_argument("--max-sites", type=int, default=40, help="most atoms in the crystal cell (default 40)")
    sc.add_argument("--limit", type=int, help="take the first N matches")
    sc.add_argument("--count", action="store_true", help="print how many crystals match and stop")
    sc.add_argument("--per-job", type=int, default=50, help="materials per job folder (default 50)")
    sc.add_argument("--checkpoint", default="~/uma/checkpoints/uma-s-1p1.pt")
    sc.add_argument("--task", default="omat")
    sc.add_argument("--adsorbates", nargs="*", default=["OH"])
    sc.add_argument("--out", type=Path, help="screen folder (default data/screens/<elements>-<count>)")
    args = parser.parse_args()
    if args.cmd == "screen":
        screen(args)
        return
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

"""Computed runs started from the app: one folder per run under data/compute/ (PE_DATA_DIR moves it), each a job folder
for the UMA worker plus meta.json. Runs wait in a queue and compute one after another (compute/runner.py), so a run
keeps going when the page closes, and a stopped run resumes from its last finished relaxation.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

from peroxide_engine.compute import jobs, surfaces

ROOT = Path(__file__).resolve().parents[3]
MODEL = {"checkpoint": "~/uma/checkpoints/uma-s-1p1.pt", "task": "omat", "name": "UMA-S 1.1, materials mode (omat)"}
ADSORBATES = ["OH", "O", "OOH", "H"]
SITES_PER_FACE = 5  # a top site per element on most oxide faces, a bridge and two hollows
SECONDS_PER_ATOM_STEP = 0.018  # median of the eight-oxide check on this computer's four performance cores
STEPS = {"slab": 30, "face": 30, "adsorbate": 70}  # median relaxation lengths in the same check, rounded up


def data_dir() -> Path:
    return Path(os.environ.get("PE_DATA_DIR", ROOT / "data")) / "compute"


def _status(folder: Path) -> dict:
    try:
        return json.loads((folder / "status.json").read_text())
    except (FileNotFoundError, json.JSONDecodeError):
        return {}


def _write_status(folder: Path, **fields) -> None:
    status = {**_status(folder), **fields, "updated": time.time()}
    tmp = folder / "status.json.tmp"
    tmp.write_text(json.dumps(status))
    os.replace(tmp, folder / "status.json")


def alive(pid: int | None) -> bool:
    if not pid:
        return False
    try:
        os.kill(pid, 0)
    except (ProcessLookupError, PermissionError):
        return False
    return True


def estimate_seconds(entries: list[dict], n_adsorbates: int) -> tuple[float, float]:
    """A rough span for the run on this computer: every face relaxed, then the adsorbates on the smallest or the
    largest face (which face wins is known only after the surface energies)."""
    total_lo = total_hi = 0.0
    for e in entries:
        faces = sum(s["n_atoms"] for s in e["slabs"]) * STEPS["slab"]
        sizes = [s["n_atoms"] * surfaces.repeats(s["cell"])[0] * surfaces.repeats(s["cell"])[1] for s in e["slabs"]]
        for size, bucket in ((min(sizes), "lo"), (max(sizes), "hi")):
            work = faces + size * STEPS["face"] + n_adsorbates * SITES_PER_FACE * (size + 3) * STEPS["adsorbate"]
            if bucket == "lo":
                total_lo += work * SECONDS_PER_ATOM_STEP
            else:
                total_hi += work * SECONDS_PER_ATOM_STEP
    return total_lo, total_hi


def prepare(rows: list[dict], processes: int | None = None) -> tuple[list[dict], list[dict]]:
    """(ready, set aside): each Materials Project row cut into its usable faces."""
    work = [(r, 240) for r in rows]
    if processes == 1 or len(rows) <= 2:
        entries = [jobs._entry_or_reason(w) for w in work]
    else:
        with ProcessPoolExecutor(processes) as pool:
            entries = list(pool.map(jobs._entry_or_reason, work))
    ready = [e for e in entries if e.get("slabs")]
    aside = [{"id": e["id"], "formula": e["formula"],
              "reason": e.get("skipped") or "every face up to Miller index 1 is polar, asymmetric or above 240 atoms"}
             for e in entries if not e.get("slabs")]
    return ready, aside


def create(name: str, ready: list[dict], aside: list[dict], adsorbates: list[str]) -> Path:
    """Write the run's folder and queue it."""
    stamp = time.strftime("%Y%m%d-%H%M%S")
    folder = data_dir() / f"{stamp}-{''.join(c if c.isalnum() else '-' for c in name)[:40].strip('-') or 'run'}"
    jobs.write_job(folder, ready, MODEL["checkpoint"], MODEL["task"], adsorbates)
    lo, hi = estimate_seconds(ready, len(adsorbates))
    meta = {"name": name, "created": time.time(), "model": MODEL["name"], "adsorbates": adsorbates,
            "materials": [{"id": e["id"], "formula": e["formula"], "label": e["label"], "faces": len(e["slabs"]),
                           **{k: e["materials_project"].get(k) for k in ("space_group", "e_hull_eV_per_atom")}}
                          for e in ready],
            "set_aside": aside, "estimate_s": [lo, hi]}
    (folder / "meta.json").write_text(json.dumps(meta, indent=1))
    _write_status(folder, state="queued", message="In the queue")
    return folder


def runs() -> list[dict]:
    """Every run, newest first, with its meta and status."""
    out = []
    root = data_dir()
    if root.exists():
        for f in sorted(root.iterdir(), reverse=True):
            if (f / "meta.json").exists():
                out.append({"folder": f, "meta": json.loads((f / "meta.json").read_text()), "status": _status(f)})
    return out


def state(run: dict) -> str:
    """queued, running, finished, stopped or failed; a worker that ended mid-run reads as stopped."""
    s = run["status"].get("state", "queued")
    if s in ("loading", "running") and not alive(run["status"].get("pid")):
        return "stopped"
    return "running" if s == "loading" else s


def runner_alive() -> bool:
    try:
        return alive(int((data_dir() / "runner.pid").read_text()))
    except (FileNotFoundError, ValueError):
        return False


def start_runner() -> None:
    """Start the queue runner unless one is already going; it outlives the page."""
    if runner_alive():
        return
    data_dir().mkdir(parents=True, exist_ok=True)
    log = (data_dir() / "runner.log").open("a")
    subprocess.Popen([sys.executable, "-m", "peroxide_engine.compute.runner"], stdout=log, stderr=subprocess.STDOUT,
                     start_new_session=True, cwd=ROOT, env={**os.environ, "PE_DATA_DIR": str(data_dir().parent)})


def stop(folder: Path) -> None:
    pid = _status(folder).get("pid")
    if alive(pid):
        os.kill(pid, 15)
    _write_status(folder, state="stopped", message="Stopped; resume picks up at the last finished relaxation")


def resume(folder: Path) -> None:
    _write_status(folder, state="queued", message="In the queue")
    start_runner()


def progress(run: dict) -> tuple[int, int]:
    """(relaxations finished, relaxations planned); the planned count assumes five sites per face until the face is
    known."""
    rows = jobs.results(run["folder"])
    done = sum(r["step"] in ("bulk", "slab", "face", "adsorbate") for r in rows)
    n_ads = len(run["meta"]["adsorbates"])
    planned = sum(2 + m["faces"] + SITES_PER_FACE * n_ads for m in run["meta"]["materials"])
    return done, max(planned, done)


def summary(folder: Path) -> list[dict]:
    """Per material: the chosen face, every face's surface energy, every site's binding energy and, with OH and
    H2O2 computed, the energy of H2O2 → 2 OH* at wide spacing."""
    rows = jobs.results(folder)
    gas = {r["molecule"]: r["energy"] for r in rows if r["step"] == "gas"}
    job = json.loads((folder / "job.json").read_text())
    out = []
    for m in job["materials"]:
        mine = [r for r in rows if r.get("material") == m["id"]]
        face = next((r for r in mine if r["step"] == "face"), None)
        sites = [r for r in mine if r["step"] == "adsorbate"]
        best = {}
        for r in sites:
            if not r["anomaly"] and (r["adsorbate"] not in best or r["E_bind"] < best[r["adsorbate"]]["E_bind"]):
                best[r["adsorbate"]] = r
        oh = best.get("OH")
        out.append({
            "id": m["id"], "formula": m["formula"], "label": m["label"],
            "state": ("error" if any(r["step"] == "error" for r in mine) else
                      "finished" if any(r["step"] == "material_done" for r in mine) else
                      "computing" if mine else "waiting"),
            "error": next((r["error"] for r in mine if r["step"] == "error"), ""),
            "face": face["hkl"] if face else None, "face_index": face["index"] if face else None,
            "face_atoms": face["n_atoms"] if face else None, "repeat": face["repeat"] if face else None,
            "faces": [{"hkl": r["hkl"], "shift": r["shift"], "atoms": r["n_atoms"], "gamma_J_m2": r["gamma_J_per_m2"],
                       "converged": r["converged"], "index": r["index"]} for r in mine if r["step"] == "slab"],
            "sites": [{"adsorbate": r["adsorbate"], "site": r["site"], "E_bind_eV": r["E_bind"],
                       "kept": not r["anomaly"], "anomaly": r["anomaly"], "converged": r["converged"],
                       "structure": f"structures/{m['id'].replace(':', '_').replace('/', '_')}__{r['adsorbate']}__"
                                    f"{r['site'].replace(' ', '-')}.cif"} for r in sites],
            "best": {a: {"site": r["site"], "E_bind_eV": r["E_bind"]} for a, r in best.items()},
            "D_eV": (2 * oh["E_minus_slab"] - gas["H2O2"]) if oh and "H2O2" in gas else None,
        })
    return out

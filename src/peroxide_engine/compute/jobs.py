"""Job folders for the UMA worker: crystals from the Materials Project, slabs cut here with pymatgen, energies computed
by worker.py in the Python that holds fairchem (PE_UMA_PYTHON), on this machine or any other.

A job folder holds everything the worker reads (job.json), so a cluster runs one folder per GPU and sends back
results.jsonl.
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

from pymatgen.core import Structure

from peroxide_engine.compute import materials, surfaces

UMA_PYTHON = os.environ.get("PE_UMA_PYTHON", "/opt/miniconda3/bin/python3")
WORKER = Path(__file__).with_name("worker.py")
SETTINGS = {"fmax": 0.05, "steps": 250, "width_A": 7.0, "nudge_A": 0.03}


def material_entry(formula: str, space_group_number: int | None = None, max_atoms: int = 240) -> dict:
    """A job entry: the Materials Project crystal (lowest above the hull, in the given space group when one is named),
    its bulk cell and every usable slab up to Miller index 1."""
    row = (materials.lowest_in_space_group(formula, space_group_number) if space_group_number
           else materials.most_stable(formula))
    return entry_from_row(row, max_atoms)


def entry_from_row(row: dict, max_atoms: int = 240) -> dict:
    st = materials.structure(row["material_id"])
    s = Structure.from_dict(st)
    usable, aside = surfaces.slabs(st, max_atoms=max_atoms)
    keep = ("hkl", "shift", "n_atoms", "area_A2", "species", "positions", "cell")
    return {"id": row["material_id"], "label": f"{row['formula']} ({row['space_group']})", "formula": row["formula"],
            "materials_project": row,
            "bulk": {"species": [site.specie.symbol for site in s], "frac": s.frac_coords.tolist(),
                     "lattice": s.lattice.matrix.tolist()},
            "slabs": [{k: r[k] for k in keep} for r in usable], "slabs_set_aside": len(aside)}


def write_job(folder: Path, entries: list[dict], checkpoint: str, task: str, adsorbates: list[str],
              gases: tuple[str, ...] = ("H2", "H2O", "H2O2")) -> Path:
    folder.mkdir(parents=True, exist_ok=True)
    job = {"checkpoint": checkpoint, "task": task, "adsorbates": adsorbates, "gases": list(gases), **SETTINGS,
           "materials": entries}
    (folder / "job.json").write_text(json.dumps(job))
    return folder


def launch(folder: Path, log: Path | None = None) -> subprocess.Popen:
    env = {**os.environ, "TORCHDYNAMO_DISABLE": "1"}
    out = (log or folder / "worker.log").open("a")
    return subprocess.Popen([UMA_PYTHON, str(WORKER), str(folder)], stdout=out, stderr=subprocess.STDOUT, env=env)


def results(folder: Path) -> list[dict]:
    path = folder / "results.jsonl"
    return [json.loads(line) for line in path.read_text().splitlines()] if path.exists() else []


def _entry_or_reason(args) -> dict:
    row, max_atoms = args
    try:
        return entry_from_row(row, max_atoms)
    except Exception as exc:  # a structure pymatgen fails to cut is listed, and the screen goes on
        return {"id": row["material_id"], "formula": row["formula"], "skipped": f"{type(exc).__name__}: {exc}"}


def prepare_screen(rows: list[dict], out: Path, checkpoint: str, task: str, adsorbates: list[str],
                   per_job: int = 50, max_atoms: int = 240, processes: int | None = None) -> dict:
    """Cut every material's faces in parallel and split the materials into job folders of `per_job` each
    (out/shard-0000, ...), with manifest.json and run_slurm.sh, the array job that runs one folder per GPU."""
    out.mkdir(parents=True, exist_ok=True)
    work = [(r, max_atoms) for r in rows]
    if processes == 1:
        entries = [_entry_or_reason(w) for w in work]
    else:
        with ProcessPoolExecutor(processes) as pool:
            entries = list(pool.map(_entry_or_reason, work, chunksize=4))
    ready = [e for e in entries if e.get("slabs")]
    aside = [{"id": e["id"], "formula": e["formula"], "reason": e.get("skipped") or "every face polar, asymmetric or "
              "above the atom limit"} for e in entries if not e.get("slabs")]
    shards = [ready[i:i + per_job] for i in range(0, len(ready), per_job)]
    for k, chunk in enumerate(shards):
        write_job(out / f"shard-{k:04d}", chunk, checkpoint, task, adsorbates)
    manifest = {"materials": len(rows), "ready": len(ready), "set_aside": aside, "shards": len(shards),
                "per_job": per_job, "checkpoint": checkpoint, "task": task, "adsorbates": adsorbates,
                "faces": sum(len(e["slabs"]) for e in ready)}
    (out / "manifest.json").write_text(json.dumps(manifest, indent=1))
    (out / "run_slurm.sh").write_text(SLURM.format(last=max(len(shards) - 1, 0)))
    shutil.copy(WORKER, out / "worker.py")  # the folder carries everything a cluster node needs besides fairchem
    return manifest


SLURM = """#!/bin/bash
# One shard per GPU: sbatch run_slurm.sh from this folder, with PE_UMA_PYTHON pointing at a Python that has fairchem
# and PE_UMA_CHECKPOINT at the UMA checkpoint (else the path in each job.json). A preempted shard resumes where it stopped when resubmitted.
#SBATCH --array=0-{last}
#SBATCH --gpus=1
#SBATCH --cpus-per-task=4
#SBATCH --time=24:00:00
SHARD=$(printf "shard-%04d" $SLURM_ARRAY_TASK_ID)
TORCHDYNAMO_DISABLE=1 ${{PE_UMA_PYTHON:-python}} worker.py "$SHARD" >> "$SHARD/worker.log" 2>&1
"""

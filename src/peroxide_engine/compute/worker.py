"""UMA worker: surface energies and adsorbate binding on prepared slabs, from one job folder, on any machine.

    TORCHDYNAMO_DISABLE=1 <python with fairchem> worker.py <job folder>

It reads <folder>/job.json (written by compute/jobs.py, which cuts the slabs with pymatgen) and imports only numpy,
ASE and fairchem, so it runs wherever UMA does: a laptop, or one GPU of a cluster per job folder. It runs on CUDA
when a GPU is present, else on the CPU. Each finished step is one JSON line in results.jsonl; a restarted job skips
every relaxation it already finished (each bulk, face and site), so a preempted cluster job picks up where it
stopped. status.json shows progress.

For each material:
1. The bulk crystal, atoms relaxed in its Materials Project cell: the energy per atom.
2. Each candidate slab (stoichiometric, nonpolar, symmetric), every atom relaxed after a 0.03 Å nudge:
   γ = (E_slab − N·E_bulk per atom) / 2A. The slab with the lowest γ is the face used.
3. That face repeated to at least `width_A` across, its lower half fixed, relaxed again: the clean reference.
4. Each adsorbate on every site the site finder returns, relaxed, set against the clean face (or the bare surface
   left under the adsorbate when that lies lower), against its gas reference:

       E_bind = E(slab + adsorbate) − E(slab) − E_ref(adsorbate)

Gas references come from the same model: OH = H2O − ½H2, O = H2O − H2, OOH = 2H2O − 3/2 H2, H = ½H2; H2O2 is
relaxed too, so H2O2 → 2OH* follows. Geometry helpers follow catalyst-engine's design/uma_worker.py.
"""

from __future__ import annotations

import hashlib
import itertools
import json
import os
import sys
import time
import traceback
from pathlib import Path

import numpy as np
from ase import Atoms
from ase.constraints import FixAtoms
from ase.geometry import find_mic
from ase.io import read, write
from ase.optimize import BFGS

NONMETALS = {"H", "B", "C", "N", "O", "F", "Si", "P", "S", "Cl", "Se", "Br", "I", "Te"}
ADSORBATE_ATOMS = {
    "H": Atoms("H"),
    "O": Atoms("O"),
    "OH": Atoms("OH", positions=[(0, 0, 0), (0, 0, 0.98)]),
    "OOH": Atoms("OOH", positions=[(0, 0, 0), (0.95, 0, 1.03), (1.86, 0, 1.40)]),
}
HEIGHTS = {  # binding atom above the site's atom height, by site kind (Å)
    "H": {"top": 1.6, "bridge": 1.2, "hollow": 1.0},
    "O": {"top": 1.9, "bridge": 1.5, "hollow": 1.2},
    "OH": {"top": 2.0, "bridge": 1.6, "hollow": 1.4},
    "OOH": {"top": 2.0, "bridge": 1.7, "hollow": 1.5},
}
_D = np.radians(112.0)  # H2O2: O-O 1.45 Å, O-H 0.97 Å, H-O-O 100°, dihedral 112°; the relaxation refines it
_c, _s = np.cos(np.radians(100.0)), np.sin(np.radians(100.0))
GAS = {
    "H2": Atoms("H2", positions=[(0, 0, 0), (0, 0, 0.74)]),
    "H2O": Atoms("OH2", positions=[(0, 0, 0), (0.76, 0, 0.59), (-0.76, 0, 0.59)]),
    "H2O2": Atoms("OOHH", positions=[(0, 0, 0), (1.45, 0, 0), (0.97 * _c, 0.97 * _s, 0),
                                     (1.45 - 0.97 * _c, 0.97 * _s * np.cos(_D), 0.97 * _s * np.sin(_D))]),
}
REFERENCE = {
    "H": {"H2": 0.5},
    "O": {"H2O": 1.0, "H2": -1.0},
    "OH": {"H2O": 1.0, "H2": -0.5},
    "OOH": {"H2O": 2.0, "H2": -1.5},
}


def relax_with(atoms: Atoms, calc, fmax: float, steps: int) -> tuple[float, int, bool]:
    atoms.calc = calc
    opt = BFGS(atoms, logfile=None)
    converged = bool(opt.run(fmax=fmax, steps=steps))
    return float(atoms.get_potential_energy()), opt.get_number_of_steps(), converged


def largest_move(before: np.ndarray, after: Atoms) -> float:
    return float(find_mic(after.positions - before, after.cell, after.pbc)[1].max())


def clean_reference(adslab: Atoms, n_slab: int, e_slab: float, energy, relax, tol: float = 0.05) -> tuple[float, float]:
    """The clean-surface energy for a binding energy: the relaxed clean slab, or, when the surface under the
    adsorbate lies lower by more than `tol` eV, that surface relaxed with the adsorbate lifted off."""
    bare = adslab[:n_slab]
    if energy(bare) >= e_slab - tol:
        return e_slab, 0.0
    e_bare = relax(bare)
    return min(e_slab, e_bare), max(0.0, e_slab - e_bare)


def _images_2d(points: np.ndarray, cell: np.ndarray) -> np.ndarray:
    shifts = [i * cell[0, :2] + j * cell[1, :2] for i in (-1, 0, 1) for j in (-1, 0, 1)]
    return np.concatenate([points + s for s in shifts])


def surface_sites(slab: Atoms, layer_A: float = 1.2, top_depth_A: float = 2.5) -> list[dict]:
    """Sites on the upper face: a top site over the highest atom of each element within `top_depth_A` of the face
    (metals first, so a cation under an oxygen layer counts), a bridge between the first outer-layer atom and its
    nearest neighbor, and up to two hollows among outer-layer atoms."""
    z = slab.positions[:, 2]
    symbols = np.array(slab.get_chemical_symbols())
    near = np.where(z > z.max() - top_depth_A)[0]
    sites: list[dict] = []
    for el in sorted(set(symbols[near]), key=lambda e: (e in NONMETALS, e)):
        i = max((i for i in near if symbols[i] == el), key=lambda i: z[i])
        sites.append({"kind": "top", "label": f"top {el}", "xy": slab.positions[i, :2].tolist(), "z": float(z[i])})
    top = np.where(z > z.max() - layer_A)[0]
    xy = slab.positions[top, :2]
    images = _images_2d(xy, slab.cell.array)
    zz = np.tile(z[top], 9)
    d = np.linalg.norm(images - xy[0], axis=1)
    d[d < 1e-6] = np.inf
    j = int(np.argmin(d))
    sites.append({"kind": "bridge", "label": "bridge", "xy": ((xy[0] + images[j]) / 2).tolist(),
                  "z": float(max(z[top[0]], zz[j]))})
    order = [k for k in np.argsort(d) if np.isfinite(d[k])][:6]
    found: list[np.ndarray] = []
    for j, k in itertools.combinations(order, 2):
        if np.linalg.norm(images[j] - images[k]) < 1.4 * min(d[j], d[k]):
            cen = (xy[0] + images[j] + images[k]) / 3
            if all(np.linalg.norm(cen - f) > 0.5 for f in found) and len(found) < 2:
                found.append(cen)
                sites.append({"kind": "hollow", "label": f"hollow {len(found)}", "xy": cen.tolist(),
                              "z": float(max(z[top[0]], zz[j], zz[k]))})
    return sites


def place(slab: Atoms, adsorbate: str, site: dict) -> Atoms:
    ads = ADSORBATE_ATOMS[adsorbate].copy()
    ads.translate([site["xy"][0], site["xy"][1], site["z"] + HEIGHTS[adsorbate][site["kind"]]])
    out = slab.copy()
    out.set_constraint(slab.constraints)
    out.extend(ads)
    return out


def anomaly(adslab: Atoms, n_slab: int, adsorbate: str) -> str:
    """'desorbed' if the adsorbate drifted off the surface, 'dissociated' if one of its bonds broke, else ''."""
    d = adslab.get_distances(n_slab, list(range(n_slab)), mic=True)
    if d.min() > 3.0:
        return "desorbed"
    ref = ADSORBATE_ATOMS[adsorbate]
    if len(ref) > 1:
        d0, d1 = ref.get_all_distances(), adslab[n_slab:].get_all_distances(mic=True)
        bonded = (d0 > 0) & (d0 < 1.5)
        if np.any(d1[bonded] > 1.5 * d0[bonded]):
            return "dissociated"
    return ""


def slab_atoms(rec: dict) -> Atoms:
    return Atoms(rec["species"], positions=rec["positions"], cell=rec["cell"], pbc=True)


def bulk_atoms(rec: dict) -> Atoms:
    return Atoms(rec["species"], scaled_positions=rec["frac"], cell=rec["lattice"], pbc=True)


def repeats(cell: np.ndarray, width_A: float) -> tuple[int, int]:
    a, b = cell[0], cell[1]
    area = np.linalg.norm(np.cross(a, b))
    return max(1, int(np.ceil(width_A * np.linalg.norm(b) / area))), max(1, int(np.ceil(width_A * np.linalg.norm(a) / area)))


def fix_lower_half(slab: Atoms) -> Atoms:
    z = slab.positions[:, 2]
    slab.set_constraint(FixAtoms(mask=z < z.min() + 0.5 * (z.max() - z.min())))
    return slab


def checkpoint_id(path: Path) -> dict:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 24), b""):
            h.update(block)
    return {"checkpoint": path.name, "sha256": h.hexdigest(), "bytes": path.stat().st_size}


def _save_status(folder: Path, **fields) -> None:
    path = folder / "status.json"
    status = json.loads(path.read_text()) if path.exists() else {}
    status.update(fields, updated=time.time())
    tmp = folder / "status.json.tmp"
    tmp.write_text(json.dumps(status))
    os.replace(tmp, path)


def _earlier(folder: Path) -> tuple[dict, dict]:
    """What an earlier run of this job finished: gas energies, and per material its bulk, slab, face and adsorbate
    rows, so a restart repeats only the relaxation that was under way."""
    gas, by_material = {}, {}
    path = folder / "results.jsonl"
    if path.exists():
        for line in path.read_text().splitlines():
            r = json.loads(line)
            if r["step"] == "gas":
                gas[r["molecule"]] = r["energy"]
            elif "material" in r and r["step"] != "error":
                m = by_material.setdefault(r["material"], {"slab": {}, "adsorbate": {}})
                if r["step"] == "slab":
                    m["slab"][r["index"]] = r
                elif r["step"] == "adsorbate":
                    m["adsorbate"][(r["adsorbate"], r["site"])] = r
                else:
                    m[r["step"]] = r
    return gas, by_material


def torch_threads() -> int:
    """CPU threads for torch: the performance cores on Apple silicon (a thread on an efficiency core holds back every
    parallel step), else all cores but two, at most 8. PE_TORCH_THREADS overrides."""
    if os.environ.get("PE_TORCH_THREADS"):
        return int(os.environ["PE_TORCH_THREADS"])
    if sys.platform == "darwin":
        try:
            return max(1, int(os.popen("sysctl -n hw.perflevel0.physicalcpu").read().strip()))
        except ValueError:
            pass
    return max(1, min(8, (os.cpu_count() or 4) - 2))


def main(folder: Path) -> None:
    job = json.loads((folder / "job.json").read_text())
    fmax, steps, width, nudge_A = job["fmax"], job["steps"], job["width_A"], job.get("nudge_A", 0.03)
    _save_status(folder, state="loading", message="Loading the UMA model", pid=os.getpid(), started=time.time())

    import torch
    from fairchem.core import FAIRChemCalculator, pretrained_mlip
    from importlib.metadata import version

    device = "cuda" if torch.cuda.is_available() else "cpu"
    if device == "cpu":
        torch.set_num_threads(torch_threads())
    ckpt = Path(os.path.expanduser(os.environ.get("PE_UMA_CHECKPOINT") or job["checkpoint"]))
    calc = FAIRChemCalculator(pretrained_mlip.load_predict_unit(str(ckpt), device=device), task_name=job["task"])
    out = (folder / "results.jsonl").open("a")

    def emit(**row) -> None:
        out.write(json.dumps(row) + "\n")
        out.flush()

    def relax(atoms: Atoms) -> tuple[float, int, bool]:
        return relax_with(atoms, calc, fmax, steps)

    def energy(atoms: Atoms) -> float:
        atoms.calc = calc
        return float(atoms.get_potential_energy())

    gas, earlier = _earlier(folder)
    emit(step="model", task=job["task"], device=device, threads=torch.get_num_threads() if device == "cpu" else None,
         fairchem=version("fairchem-core"), torch=torch.__version__, **checkpoint_id(ckpt))
    _save_status(folder, state="running", message="Gas-phase molecules")
    for name in job["gases"]:
        if name in gas:
            continue
        mol = GAS[name].copy()
        mol.cell = [15.0, 15.1, 15.2]
        mol.center()
        mol.pbc = True
        e, n, conv = relax(mol)
        gas[name] = e
        emit(step="gas", molecule=name, energy=e, steps=n, converged=conv)
    refs = {a: sum(w * gas[g] for g, w in REFERENCE[a].items()) for a in job["adsorbates"]}
    (folder / "structures").mkdir(exist_ok=True)

    for k, m in enumerate(job["materials"]):
        before = earlier.get(m["id"], {"slab": {}, "adsorbate": {}})
        if "material_done" in before:
            continue
        tag, label = m["id"].replace(":", "_").replace("/", "_"), m["label"]
        try:
            _save_status(folder, message=f"{label}: bulk", material=k + 1, materials=len(job["materials"]))
            if "bulk" in before:
                per_atom = before["bulk"]["per_atom"]
            else:
                bulk = bulk_atoms(m["bulk"])
                bulk.rattle(stdev=nudge_A, seed=0)
                e_bulk, n, conv = relax(bulk)
                per_atom = e_bulk / len(bulk)
                emit(step="bulk", material=m["id"], energy=e_bulk, n_atoms=len(bulk), per_atom=per_atom, steps=n,
                     converged=conv)
            faces = []
            for i, rec in enumerate(m["slabs"]):
                if i in before["slab"]:
                    faces.append((before["slab"][i]["gamma_eV_per_A2"], i))
                    continue
                _save_status(folder, message=f"{label}: face {rec['hkl']} ({i + 1} of {len(m['slabs'])})")
                t0 = time.time()
                s = slab_atoms(rec)
                built = s.positions.copy()
                s.rattle(stdev=nudge_A, seed=0)
                e, n, conv = relax(s)
                area = float(np.linalg.norm(np.cross(s.cell.array[0], s.cell.array[1])))
                gamma = (e - len(s) * per_atom) / (2 * area)
                write(folder / "structures" / f"{tag}__slab{i}.cif", s)
                faces.append((gamma, i))
                emit(step="slab", material=m["id"], index=i, hkl=rec["hkl"], shift=rec["shift"], n_atoms=len(s),
                     area_A2=area, energy=e, gamma_eV_per_A2=gamma, gamma_J_per_m2=gamma * 16.0218,
                     moved_A=round(largest_move(built, s), 3), steps=n, converged=conv,
                     seconds=round(time.time() - t0, 1))
            gamma, i = min(faces)
            best = read(folder / "structures" / f"{tag}__slab{i}.cif")
            ra, rb = repeats(best.cell.array, width)
            slab = fix_lower_half(best.repeat((ra, rb, 1)))
            if "face" in before:
                slab.positions = read(folder / "structures" / f"{tag}__face.cif").positions
                e_slab = before["face"]["energy"]
            else:
                _save_status(folder, message=f"{label}: clean face {m['slabs'][i]['hkl']} at {ra}×{rb}")
                e_slab, n, conv = relax(slab)
                write(folder / "structures" / f"{tag}__face.cif", slab)
                emit(step="face", material=m["id"], index=i, hkl=m["slabs"][i]["hkl"], repeat=[ra, rb],
                     n_atoms=len(slab), energy=e_slab, steps=n, converged=conv)
            sites = surface_sites(slab)
            for a in job["adsorbates"]:
                for site in sites:
                    if (a, site["label"]) in before["adsorbate"]:
                        continue
                    _save_status(folder, message=f"{label}: {a} on the {site['label']} site")
                    t0 = time.time()
                    adslab = place(slab, a, site)
                    adslab.rattle(stdev=nudge_A, seed=1)
                    e, n, conv = relax(adslab)
                    e_ref, drop = clean_reference(adslab, len(slab), e_slab, energy, lambda x: relax(x)[0])
                    name = f"{tag}__{a}__{site['label'].replace(' ', '-')}.cif"
                    write(folder / "structures" / name, adslab)
                    emit(step="adsorbate", material=m["id"], adsorbate=a, site=site["label"], energy=e,
                         E_bind=e - e_ref - refs[a], E_minus_slab=e - e_ref, slab_drop_eV=round(drop, 3),
                         anomaly=anomaly(adslab, len(slab), a), steps=n, converged=conv,
                         seconds=round(time.time() - t0, 1))
            emit(step="material_done", material=m["id"])
        except Exception as exc:  # one material failing leaves the rest of the job running
            emit(step="error", material=m["id"], error=f"{type(exc).__name__}: {exc}")
            traceback.print_exc()
    out.close()
    _save_status(folder, state="finished", message="Finished", finished=time.time())


if __name__ == "__main__":
    root = Path(sys.argv[1])
    try:
        main(root)
    except Exception as exc:
        traceback.print_exc()
        _save_status(root, state="failed", message=f"{type(exc).__name__}: {exc}")
        raise

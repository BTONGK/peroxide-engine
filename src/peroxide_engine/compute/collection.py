"""Every computed result in one place: runs started in the app (data/compute/), screens prepared with `pe screen`
(data/screens/*/shard-*), and the model checks' runs (experiments/*/runs/*). Each row names the model and the run
that made it, so the same crystal computed twice appears twice, once per run.
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from peroxide_engine.compute import screens

ROOT = screens.ROOT
MODELS = {("uma-s-1p1.pt", "omat"): "UMA-S 1.1, materials mode", ("uma-s-1p2p1.pt", "oc22"): "UMA-S 1.2.1, oxide mode"}


def sources() -> list[dict]:
    """Every job folder holding results, with where it came from."""
    out = []
    for run in screens.runs():
        out.append({"folder": run["folder"], "run": run["meta"]["name"], "origin": "Compute a catalyst",
                    "created": run["meta"]["created"]})
    data = screens.data_dir().parent
    for shard in sorted((data / "screens").glob("*/shard-*")):
        out.append({"folder": shard, "run": f"{shard.parent.name} / {shard.name}", "origin": "Screen (pe screen)",
                    "created": (shard / "job.json").stat().st_mtime})
    for folder in sorted((ROOT / "experiments").glob("*/runs/*")):
        out.append({"folder": folder, "run": f"{folder.parent.parent.name} / {folder.name}", "origin": "Model check",
                    "created": (folder / "job.json").stat().st_mtime})
    return [s for s in out if (s["folder"] / "results.jsonl").exists() and (s["folder"] / "job.json").exists()]


def fingerprint() -> tuple:
    """Changes whenever any results file grows, for caching the tables."""
    return tuple((str(s["folder"]), (s["folder"] / "results.jsonl").stat().st_mtime) for s in sources())


def model_name(job: dict) -> str:
    ckpt = Path(job["checkpoint"]).name
    return MODELS.get((ckpt, job["task"]), f"{ckpt} ({job['task']})")


def tables() -> dict[str, pd.DataFrame]:
    """crystals (one row per crystal per run), faces (every face relaxed) and sites (every binding energy)."""
    crystals, faces, sites = [], [], []
    for src in sources():
        job = json.loads((src["folder"] / "job.json").read_text())
        mp = {m["id"]: m.get("materials_project", {}) for m in job["materials"]}
        common = {"Model": model_name(job), "Run": src["run"], "Origin": src["origin"]}
        for m in screens.summary(src["folder"]):
            if not m["faces"]:
                continue
            row = mp.get(m["id"], {})
            gamma = min(f["gamma_J_m2"] for f in m["faces"])
            crystal = {"Formula": m["formula"], "Materials Project": m["id"], "Space group": row.get("space_group"),
                       "Above the hull (eV/atom)": row.get("e_hull_eV_per_atom"), "State": m["state"],
                       "Face": m["face"] and "(" + " ".join(map(str, m["face"])) + ")",
                       "Surface energy (J/m²)": gamma, "Faces computed": len(m["faces"])}
            for a in job["adsorbates"]:
                b = m["best"].get(a)
                crystal[f"{a} binding (eV)"] = b and b["E_bind_eV"]
                crystal[f"{a} site"] = b and b["site"]
            if "OH" in job["adsorbates"]:
                crystal["H2O2 → 2 OH* (eV)"] = m["D_eV"]
            crystals.append({**crystal, **common})
            for f in m["faces"]:
                faces.append({"Formula": m["formula"], "Materials Project": m["id"],
                              "Face": "(" + " ".join(map(str, f["hkl"])) + ")", "Termination shift": f["shift"],
                              "Atoms": f["atoms"], "Surface energy (J/m²)": f["gamma_J_m2"],
                              "Lowest of its crystal": f["index"] == m["face_index"],
                              "Converged": f["converged"], **common})
            for s in m["sites"]:
                sites.append({"Formula": m["formula"], "Materials Project": m["id"], "Face": crystal["Face"],
                              "Molecule": s["adsorbate"], "Site": s["site"], "Binding (eV)": s["E_bind_eV"],
                              "Kept": s["kept"], "Flag": s["anomaly"], "Converged": s["converged"], **common})
    return {"crystals": pd.DataFrame(crystals), "faces": pd.DataFrame(faces), "sites": pd.DataFrame(sites)}


def measured() -> pd.DataFrame:
    """The measured values the computed data are tested against: the oxide check's test set, then the library's
    ready catalysts."""
    import yaml

    from peroxide_engine.catalysts import library

    rows = []
    check = yaml.safe_load((ROOT / "experiments" / "oxide_descriptor_check" / "data.yaml").read_text())
    for ox in check["oxides"]:
        rows.append({"Catalyst": f"{ox['formula']} ({ox['phase']})", "Quantity": "Activation energy",
                     "Value": ox["Ea_kJ"], "Units": "kJ/mol", "Source": "Lousada 2013, Table 4",
                     "Used in": "Oxide check (test set)"})
        if "k2_m_s" in ox:
            rows.append({"Catalyst": f"{ox['formula']} ({ox['phase']})", "Quantity": "Rate per surface area k2",
                         "Value": ox["k2_m_s"], "Units": "m/s",
                         "Source": "Lousada 2013, Table 3" + (f" ({ox['k2_note']})" if "k2_note" in ox else ""),
                         "Used in": "Oxide check (test set)"})
    for e in library.runnable():
        law = e.law
        for label, v in (("Rate per surface area k2", law.k2), ("Intercept", law.intercept), ("Activation energy", law.Ea),
                         ("Specific surface area", law.specific_area)):
            if v is not None:
                rows.append({"Catalyst": e.name, "Quantity": label, "Value": v.printed, "Units": v.units,
                             "Source": f"{e.source.cite}, {v.locator}", "Used in": "Catalyst library (Run page)"})
    return pd.DataFrame(rows)

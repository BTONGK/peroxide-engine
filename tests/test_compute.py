"""The computed-catalyst pipeline: surfaces cut from a crystal, the worker's geometry, the check's ranking statistics."""

import importlib.util
from pathlib import Path

import numpy as np
from pymatgen.core import Lattice, Structure

from peroxide_engine.compute import surfaces, worker

ROOT = Path(__file__).resolve().parents[1]


def ceria() -> dict:
    """Fluorite CeO2, a = 5.41 Å, as a pymatgen dict (the form the Materials Project serves)."""
    return Structure.from_spacegroup("Fm-3m", Lattice.cubic(5.41), ["Ce", "O"], [[0, 0, 0], [0.25, 0.25, 0.25]]).as_dict()


def test_usable_slabs_are_stoichiometric_nonpolar_and_symmetric():
    usable, aside = surfaces.slabs(ceria())
    assert {tuple(s["hkl"]) for s in usable} >= {(1, 1, 1)}
    assert all(s["stoichiometric"] and s["symmetric"] and not s["polar"] for s in usable)
    assert all(not (s["stoichiometric"] and s["symmetric"] and not s["polar"]) or s["n_atoms"] > 240 for s in aside)


def test_repeats_reach_the_width_across_a_skewed_cell():
    hexagonal = [[3.83, 0, 0], [3.83 * 0.5, 3.83 * np.sqrt(3) / 2, 0], [0, 0, 30]]
    assert surfaces.repeats(hexagonal, 7.0) == worker.repeats(np.array(hexagonal), 7.0) == (3, 3)


def test_site_finder_reaches_a_metal_under_an_oxygen_layer():
    usable, _ = surfaces.slabs(ceria())
    rec = next(s for s in usable if tuple(s["hkl"]) == (1, 1, 1))
    slab = worker.slab_atoms(rec)
    z = slab.positions[:, 2]
    assert slab[int(np.argmax(z))].symbol == "O"  # CeO2(111) ends in oxygen
    labels = [s["label"] for s in worker.surface_sites(slab)]
    assert "top Ce" in labels and "top O" in labels and "bridge" in labels


def test_gas_geometry_and_anomaly_flags():
    h2o2 = worker.GAS["H2O2"]
    assert np.isclose(h2o2.get_distance(0, 1), 1.45) and np.isclose(h2o2.get_distance(0, 2), 0.97)
    assert np.isclose(h2o2.get_distance(1, 3), 0.97, atol=1e-6)
    usable, _ = surfaces.slabs(ceria())
    slab = worker.slab_atoms(next(s for s in usable if tuple(s["hkl"]) == (1, 1, 1)))
    site = worker.surface_sites(slab)[0]
    whole = worker.place(slab, "OH", site)
    assert worker.anomaly(whole, len(slab), "OH") == ""
    broken = whole.copy()
    broken.positions[-1, 2] += 1.0
    assert worker.anomaly(broken, len(slab), "OH") == "dissociated"


def test_ranking_statistics():
    spec = importlib.util.spec_from_file_location("odc", ROOT / "experiments" / "oxide_descriptor_check" / "check.py")
    odc = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(odc)
    x = [1, 2, 3, 4, 5]
    assert np.isclose(odc.spearman(x, [10, 20, 30, 40, 50]), 1.0)
    assert np.isclose(odc.exact_p(x, [10, 20, 30, 40, 50], +1), 1 / 120)
    assert np.isclose(odc.exact_p(x, [50, 40, 30, 20, 10], -1), 1 / 120)
    assert np.allclose(odc.loo([0, 1, 2, 3], [1, 3, 5, 7]), 0.0, atol=1e-9)


def test_search_reads_every_page_and_caches_structures(monkeypatch, tmp_path):
    from peroxide_engine.compute import materials

    pages = [[{"material_id": f"mp-{i}", "formula_pretty": "CeO2", "structure": ceria()} for i in range(k, k + n)]
             for k, n in ((0, 3), (3, 3), (6, 1))]
    seen = []

    def fake_get(params):
        seen.append(params["_skip"])
        return pages[params["_skip"] // 3]

    monkeypatch.setattr(materials, "_get", fake_get)
    monkeypatch.setattr(materials, "CACHE", tmp_path)
    rows = materials.search(elements=["O"], per_page=3)
    assert [r["material_id"] for r in rows] == [f"mp-{i}" for i in range(7)] and seen == [0, 3, 6]
    assert len(list(tmp_path.glob("mp-*.json"))) == 7
    assert len(materials.search(elements=["O"], per_page=3, limit=4)) == 4


def test_prepare_screen_writes_self_contained_shards(monkeypatch, tmp_path):
    from peroxide_engine.compute import jobs, materials

    monkeypatch.setattr(materials, "structure", lambda mid: ceria())
    rows = [{"material_id": f"mp-{i}", "formula": "CeO2", "space_group": "Fm-3m"} for i in range(3)]
    m = jobs.prepare_screen(rows, tmp_path, "~/uma/checkpoints/uma-s-1p1.pt", "omat", ["OH"], per_job=2, processes=1)
    assert m["ready"] == 3 and m["shards"] == 2
    assert (tmp_path / "shard-0000" / "job.json").exists() and (tmp_path / "shard-0001" / "job.json").exists()
    assert (tmp_path / "worker.py").read_text() == (ROOT / "src" / "peroxide_engine" / "compute" / "worker.py").read_text()
    assert "#SBATCH --array=0-1" in (tmp_path / "run_slurm.sh").read_text()


def test_a_run_queues_reports_progress_and_joins_the_collection(monkeypatch, tmp_path):
    import json

    from peroxide_engine.compute import collection, materials, screens

    monkeypatch.setenv("PE_DATA_DIR", str(tmp_path))
    monkeypatch.setattr(materials, "structure", lambda mid: ceria())
    row = {"material_id": "mp-ce", "formula": "CeO2", "space_group": "Fm-3m", "e_hull_eV_per_atom": 0.0}
    ready, aside = screens.prepare([row], processes=1)
    folder = screens.create("Ceria test", ready, aside, ["OH"])
    run = screens.runs()[0]
    assert screens.state(run) == "queued" and screens.progress(run)[0] == 0
    assert screens.estimate_seconds(ready, 1)[0] > 0
    lines = [{"step": "gas", "molecule": "H2O2", "energy": -17.0},
             {"step": "bulk", "material": "mp-ce", "per_atom": -8.0},
             {"step": "slab", "material": "mp-ce", "index": 0, "hkl": [1, 1, 1], "shift": 0.1, "n_atoms": 18,
              "gamma_J_per_m2": 0.67, "converged": True},
             {"step": "face", "material": "mp-ce", "index": 0, "hkl": [1, 1, 1], "repeat": [3, 3], "n_atoms": 162,
              "energy": -1400.0},
             {"step": "adsorbate", "material": "mp-ce", "adsorbate": "OH", "site": "top Ce", "E_bind": 1.6,
              "E_minus_slab": -9.0, "anomaly": "", "converged": True},
             {"step": "adsorbate", "material": "mp-ce", "adsorbate": "OH", "site": "top O", "E_bind": 1.2,
              "E_minus_slab": -9.4, "anomaly": "dissociated", "converged": True},
             {"step": "material_done", "material": "mp-ce"}]
    (folder / "results.jsonl").write_text("".join(json.dumps(r) + "\n" for r in lines))
    m = screens.summary(folder)[0]
    assert m["state"] == "finished" and m["best"]["OH"]["site"] == "top Ce"  # the dissociated site stays out
    assert np.isclose(m["D_eV"], 2 * -9.0 + 17.0)
    t = collection.tables()
    ours = t["crystals"][t["crystals"]["Origin"] == "Compute a catalyst"]
    assert list(ours["Formula"]) == ["CeO2"] and len(t["sites"][t["sites"]["Run"] == "Ceria test"]) == 2

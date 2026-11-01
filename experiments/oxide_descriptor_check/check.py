"""Computed OH binding on eight oxides against their measured activity, under PASS_RULE.md. Run:
uv run python experiments/oxide_descriptor_check/prepare.py     (writes the UMA job folders in runs/)
TORCHDYNAMO_DISABLE=1 /opt/miniconda3/bin/python3 src/peroxide_engine/compute/worker.py experiments/oxide_descriptor_check/runs/<model>
uv run python experiments/oxide_descriptor_check/check.py
"""

import itertools
import json
import math
from pathlib import Path

import numpy as np
import yaml

from peroxide_engine import checks
from peroxide_engine.compute import jobs

HERE = Path(__file__).parent
data = yaml.safe_load((HERE / "data.yaml").read_text())
PRIMARY, SECONDARY = "uma-s-1p1_omat", "uma-s-1p2p1_oc22"


def number(printed: str) -> float:
    return float(str(printed).split("±")[0])


def ranks(x) -> np.ndarray:
    x = np.asarray(x, float)
    order = np.argsort(x, kind="stable")
    r = np.empty(len(x))
    r[order] = np.arange(1, len(x) + 1)
    for v in np.unique(x):  # ties share their mean rank
        r[x == v] = r[x == v].mean()
    return r


def spearman(x, y) -> float:
    a, b = ranks(x), ranks(y)
    return float(np.corrcoef(a, b)[0, 1])


def exact_p(x, y, sign: int) -> float:
    """Share of all orderings of y whose rho with x is at least as far in the stated direction (one-sided)."""
    rho = sign * spearman(x, y)
    a, b = ranks(x), ranks(y)
    hits = total = 0
    for perm in itertools.permutations(b):
        total += 1
        hits += sign * np.corrcoef(a, perm)[0, 1] >= rho - 1e-12
    return hits / total


def descriptor(rows: list[dict]) -> dict:
    """Per oxide: the face, its surface energy, every OH site, and D = 2[E(slab + OH) − E(slab)] − E(H2O2)."""
    gas = {r["molecule"]: r["energy"] for r in rows if r["step"] == "gas"}
    out = {}
    for r in rows:
        if r["step"] == "face":
            out[r["material"]] = dict(hkl=r["hkl"], repeat=r["repeat"], n_atoms=r["n_atoms"], sites=[])
    for r in rows:
        if r["step"] == "slab" and r["material"] in out:
            out[r["material"]].setdefault("gammas", []).append((r["hkl"], r["gamma_J_per_m2"]))
        if r["step"] == "adsorbate" and r["adsorbate"] == "OH":
            out[r["material"]]["sites"].append(dict(site=r["site"], D=2 * r["E_minus_slab"] - gas["H2O2"],
                                                    E_bind_OH=r["E_bind"], anomaly=r["anomaly"],
                                                    converged=r["converged"]))
    for v in out.values():
        whole = [s for s in v["sites"] if not s["anomaly"]]
        best = min(whole, key=lambda s: s["D"]) if whole else None
        v.update(D=best["D"] if best else math.nan, site=best["site"] if best else "",
                 gamma=min(g for _, g in v["gammas"]), whole_sites=len(whole))
    return out


def loo(D, Ea) -> list[float]:
    D, Ea = np.asarray(D), np.asarray(Ea)
    errors = []
    for i in range(len(D)):
        keep = np.arange(len(D)) != i
        a, b = np.polyfit(D[keep], Ea[keep], 1)
        errors.append(float(a * D[i] + b - Ea[i]))
    return errors


def score(model: str) -> dict:
    rows = jobs.results(HERE / "runs" / model)
    job = json.loads((HERE / "runs" / model / "job.json").read_text())
    by_formula = {m["formula"]: m["id"] for m in job["materials"]}
    desc = descriptor(rows)
    table = []
    for ox in data["oxides"]:
        mid = by_formula[ox["formula"]]
        d = desc.get(mid, {})
        table.append(dict(oxide=ox["formula"], phase=ox["phase"], material_id=mid, face=str(d.get("hkl", "")),
                          surface_energy_J_m2=d.get("gamma", math.nan), site=d.get("site", ""),
                          D_eV=d.get("D", math.nan), whole_sites=d.get("whole_sites", 0),
                          Ea_measured_kJ=number(ox["Ea_kJ"]),
                          k2_measured_m_s=number(ox["k2_m_s"]) if "k2_m_s" in ox else math.nan))
    done = [t for t in table if t["D_eV"] == t["D_eV"]]
    out = dict(model=model, n=len(done), rows=table)
    if len(done) < len(table):
        return out
    D = [t["D_eV"] for t in table]
    Ea = [t["Ea_measured_kJ"] for t in table]
    rate = [(t["D_eV"], math.log10(t["k2_measured_m_s"])) for t in table if t["k2_measured_m_s"] == t["k2_measured_m_s"]]
    errs = loo(D, Ea)
    out.update(
        rho_Ea=spearman(D, Ea), p_Ea=exact_p(D, Ea, +1), p_Ea_opposite=exact_p(D, Ea, -1),
        rho_rate=spearman(*zip(*rate)), p_rate=exact_p(*zip(*rate), -1), n_rate=len(rate),
        loo_within_6=int(sum(abs(e) <= 6 for e in errs)), loo_median_abs_kJ=float(np.median(np.abs(errs))),
        loo_errors_kJ=errs, slope_kJ_per_eV=float(np.polyfit(D, Ea, 1)[0]))
    out.update(passes_Ea=bool(out["n"] >= 5 and out["rho_Ea"] >= 0.6 and out["p_Ea"] <= 0.05),
               passes_rate=bool(out["n_rate"] >= 5 and out["rho_rate"] <= -0.6 and out["p_rate"] <= 0.05),
               passes_value=bool(out["loo_within_6"] >= 6))
    return out


def main():
    primary, secondary = score(PRIMARY), score(SECONDARY)
    if "passes_Ea" not in primary:
        print(f"{PRIMARY}: {primary['n']} of {len(primary['rows'])} oxides computed; the summary waits for all")
        return
    verdict = lambda ok: "passed" if ok else "missed"
    summary = dict(
        verdict=verdict(primary["passes_Ea"]), rate_verdict=verdict(primary["passes_rate"]),
        value_verdict=verdict(primary["passes_value"]),
        n=primary["n"], rho_Ea=primary["rho_Ea"], p_Ea=primary["p_Ea"], p_Ea_opposite=primary["p_Ea_opposite"],
        n_rate=primary["n_rate"],
        rho_rate=primary["rho_rate"], p_rate=primary["p_rate"], loo_within_6=primary["loo_within_6"],
        loo_median_abs_kJ=primary["loo_median_abs_kJ"], slope_kJ_per_eV=primary["slope_kJ_per_eV"],
        secondary_model=SECONDARY, secondary_state="computed" if "rho_Ea" in secondary else "held",
        secondary_rho_Ea=secondary.get("rho_Ea", math.nan), secondary_rho_rate=secondary.get("rho_rate", math.nan),
        rows=primary["rows"], **({"secondary_rows": secondary["rows"]} if "rho_Ea" in secondary else {}),
    )
    checks.write_summary(HERE, summary)
    for t in primary["rows"]:
        print(f"{t['oxide']:<6} face {t['face']:<11} γ {t['surface_energy_J_m2']:.2f} J/m²  D {t['D_eV']:+.3f} eV"
              f"  ({t['site']}, {t['whole_sites']} whole)  Ea {t['Ea_measured_kJ']:g}")
    print(f"Ea ranking rho {primary['rho_Ea']:+.3f} p {primary['p_Ea']:.4f} -> {summary['verdict']}")
    print(f"rate ranking rho {primary['rho_rate']:+.3f} p {primary['p_rate']:.4f} -> {summary['rate_verdict']}")
    print(f"Ea value: {primary['loo_within_6']} of 8 within 6 kJ/mol -> {summary['value_verdict']}")
    if "rho_Ea" in secondary:
        print(f"{SECONDARY}: Ea rho {secondary['rho_Ea']:+.3f} p {secondary['p_Ea']:.4f}; "
              f"rate rho {secondary['rho_rate']:+.3f}")


if __name__ == "__main__":
    main()

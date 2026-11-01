"""Heat released by H2O2 solutions: NBS dilution series against Solvay's published table. Run: uv run python
experiments/decomposition_heat/check.py"""

from pathlib import Path

import numpy as np
import yaml
from scipy.interpolate import PchipInterpolator

from peroxide_engine import checks, sources, thermo

HERE = Path(__file__).parent
data = yaml.safe_load((HERE / "data.yaml").read_text())
btu = sources.value("btu_it_per_pound") / 1000  # kJ/kg per Btu/lb
kcal = sources.value("thermochemical_calorie")  # kJ per kcal
M = thermo.molar_mass_h2o2()


def half_unit(printed: str) -> float:
    decimals = len(printed.split(".")[1]) if "." in printed else 0
    return 0.5 * 10.0**-decimals


rows = []
for item in data["per_pound_of_solution"]:
    w, printed = item["mass_fraction"], item["printed"]
    mol_per_kg = 1000 / M if w == 1.0 else w * 1000 / M
    q = thermo.decomposition_heat_of_solution(w).value  # kJ per kg of solution (or of H2O2 at w = 1)
    computed = q / btu
    allowed = half_unit(printed) + 0.8 * mol_per_kg / btu
    rows.append(dict(kind="per pound", mass_fraction=w, printed=float(printed), computed=computed,
                     difference=computed - float(printed), allowed=allowed,
                     within=bool(abs(computed - float(printed)) <= allowed), units="Btu/lb"))
for item in data["per_mole"]:
    w, printed = item["mass_fraction"], item["printed"]
    q = thermo.decomposition_heat_of_solution(w).value / (w * 1000 / M) / kcal
    allowed = half_unit(printed) + 0.8 / kcal
    rows.append(dict(kind="per mole", mass_fraction=w, printed=float(printed), computed=q,
                     difference=q - float(printed), allowed=allowed,
                     within=bool(abs(q - float(printed)) <= allowed), units="kcal/mol"))

# leave-one-out interpolation error of the dilution series
u, h = thermo._U, thermo._H
loo = []
for i in range(1, len(u) - 1):
    keep = np.arange(len(u)) != i
    pred = float(PchipInterpolator(u[keep], h[keep])(u[i]))
    loo.append(dict(n_water=float(u[i] / (1 - u[i])), tabulated=float(h[i]), predicted=pred, error=pred - float(h[i])))
loo_max = max(abs(x["error"]) for x in loo)
# the error that matters is between neighbours; the leave-one-out gap is twice the table spacing, so it bounds it
per_pound = [r for r in rows if r["kind"] == "per pound"]
passed = all(r["within"] for r in per_pound) and loo_max < 0.1
summary = dict(
    verdict="passed" if passed else "missed",
    per_pound_within=sum(r["within"] for r in per_pound), per_pound_total=len(per_pound),
    per_mole_within=sum(r["within"] for r in rows if r["kind"] == "per mole"),
    per_mole_total=sum(r["kind"] == "per mole" for r in rows),
    loo_max_kJ=loo_max,
    dilute_dH=thermo.decomposition("aq")["dH"].value, pure_dH=thermo.decomposition("l")["dH"].value,
    rows=rows, leave_one_out=loo,
)
print(checks.write_summary(HERE, summary))
for r in rows:
    print(f"{r['kind']:9s} {r['mass_fraction']:.2f} printed {r['printed']:8.2f} computed {r['computed']:8.2f} "
          f"diff {r['difference']:+.2f} allowed {r['allowed']:.2f} {r['units']} {'ok' if r['within'] else 'MISS'}")
print("leave-one-out max error", round(loo_max, 4), "kJ/mol;", summary["verdict"])

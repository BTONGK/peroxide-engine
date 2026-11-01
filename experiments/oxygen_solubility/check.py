"""Oxygen solubility: Benson-Krause-Peterson against Battino's recommended values and the USGS air-saturation formula.
Run: uv run python experiments/oxygen_solubility/check.py"""

from pathlib import Path

import numpy as np
import yaml

from peroxide_engine import checks, oxygen, sources, thermo

HERE = Path(__file__).parent
data = yaml.safe_load((HERE / "data.yaml").read_text())
T0 = sources.value("celsius_zero")
ATM = sources.value("standard_atmosphere")

rows = []
for T, x_rec in data["recommended_1e5_x"].items():
    x = float(oxygen.mole_fraction_1atm(T)) * 1e5
    in_range = 278.15 <= T <= 328.15
    diff = x / x_rec - 1
    rows.append(dict(T_K=T, t_C=T - T0, recommended=x_rec, computed=x, rel_diff=diff, in_app_range=in_range,
                     within=bool(abs(diff) <= 0.005) if in_range else None))

c = sources.coefficients("o2_air_saturation_bk1984")
air = []
for t_C in [0.01, 5, 10, 15, 20, 25, 30, 35, 40]:  # IAPWS water starts at the triple point
    T = T0 + t_C
    usgs = float(np.exp(c["c0"] + c["c1"] / T + c["c2"] / T**2 + c["c3"] / T**3 + c["c4"] / T**4))
    mine = float(oxygen.air_saturation(T, ATM)) * oxygen.molar_mass() * 1000
    air.append(dict(t_C=t_C, usgs_mg_L=usgs, computed_mg_L=mine, rel_diff=mine / usgs - 1,
                    within=bool(abs(mine / usgs - 1) <= 0.001)))

dH = float(oxygen.solution_enthalpy(298.15))
in_range = [r for r in rows if r["in_app_range"]]
passed = all(r["within"] for r in in_range) and all(a["within"] for a in air)
summary = dict(
    verdict="passed" if passed else "missed",
    max_rel_diff_app_range=max(abs(r["rel_diff"]) for r in in_range),
    max_rel_diff_above=max(abs(r["rel_diff"]) for r in rows if r["T_K"] > 328.15),
    max_rel_diff_air=max(abs(a["rel_diff"]) for a in air),
    air_sat_25C_mg_L=next(a["computed_mg_L"] for a in air if a["t_C"] == 25),
    dH_solution_kJ=dH, dH_battino_kJ=data["enthalpy_of_solution_298_kJ"],
    dH_nbs_kJ=thermo.nbs("O2(ao)", "dfH"),
    rows=rows, air=air,
)
print(checks.write_summary(HERE, summary))
for r in rows:
    print(f"{r['t_C']:5.1f} °C  rec {r['recommended']:.3f}  BKP {r['computed']:.4f}  {r['rel_diff']*100:+.2f}%  {r['within']}")
for a in air:
    print(f"{a['t_C']:5.2f} °C  USGS {a['usgs_mg_L']:.3f}  engine {a['computed_mg_L']:.3f} mg/L  {a['rel_diff']*100:+.3f}%  {a['within']}")
print("dH solution", round(dH, 3), "kJ/mol;", summary["verdict"])

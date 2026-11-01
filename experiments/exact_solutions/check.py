"""Exact solutions and conservation for the vessel model and the fits. Run: uv run python
experiments/exact_solutions/check.py"""

from pathlib import Path

import numpy as np

from peroxide_engine import analysis, checks, oxygen, thermo, vessel, water

HERE = Path(__file__).parent
k, c0 = 2e-3, 0.2
results = {}

# 1 and 2
worst_exact, worst_balance = 0.0, 0.0
modes = []
for hs in ("burette", "closed", "open"):
    for kla in (None, 0.05):
        for custom in (False, True):
            if custom and kla is None and hs != "burette":
                continue
            s = vessel.Setup(c0_M=c0, k_per_s=k, headspace=hs, kLa_per_s=kla)
            run = vessel.simulate(s, rate=(lambda c: k * c) if custom else None)
            t = run.series["time_s"].to_numpy()
            exact = c0 * np.exp(-k * t) * 1e3
            err = float(np.max(np.abs(run.series["H2O2_mM"].to_numpy() / exact - 1)))
            b = run.balance
            total = b["start"] + b["formed"]
            imb = float(np.max(np.abs((b["liquid"] + b["gas"] + b["released"] - total) / total)))
            worst_exact, worst_balance = max(worst_exact, err), max(worst_balance, imb)
            modes.append(dict(headspace=hs, kLa=kla, rate_function=custom, exact_rel_err=err, balance_rel_err=imb))
results["exact_ok"] = worst_exact <= 1e-7
results["balance_ok"] = worst_balance <= 1e-9

# 3
conv = []
for hs, col in (("burette", "gas_collected_mL"), ("closed", "pressure_kPa")):
    eq = vessel.simulate(vessel.Setup(c0_M=c0, k_per_s=k, headspace=hs))
    fast = vessel.simulate(vessel.Setup(c0_M=c0, k_per_s=k, headspace=hs, kLa_per_s=1000 * k))
    a, b_ = eq.series[col].iloc[-1], fast.series[col].iloc[-1]
    conv.append(dict(headspace=hs, readout=col, equilibrium=a, fast_kLa=b_, rel_diff=b_ / a - 1))
results["limit_ok"] = all(abs(x["rel_diff"]) <= 1e-3 for x in conv)

# 4
heat = []
for hs in ("burette", "closed"):
    for kla in (None, 0.05):
        run = vessel.simulate(vessel.Setup(c0_M=c0, k_per_s=k, headspace=hs, kLa_per_s=kla))
        b, VL = run.balance, run.setup.liquid_mL / 1000
        decomposed = (c0 - run.series["H2O2_mM"].iloc[-1] / 1e3) * VL
        moved = (b["gas"].iloc[-1] - b["gas"].iloc[0])
        T, P = run.setup.T_K, run.setup.pressure_Pa
        pw = float(water.vapor_pressure(T))
        carried = pw / (P - pw) if hs == "burette" else 0.0
        dvap = thermo.nbs("H2O(g)", "dfH") - thermo.nbs("H2O(l)", "dfH")
        expect = (-thermo.decomposition_enthalpy_to_dissolved_o2(T) * decomposed
                  + float(oxygen.solution_enthalpy(T)) * moved - carried * dvap * moved) * 1e3
        got = run.series["heat_J"].iloc[-1]
        heat.append(dict(headspace=hs, kLa=kla, bookkeeping_J=expect, computed_J=got, rel_diff=got / expect - 1))
results["heat_ok"] = all(abs(h["rel_diff"]) <= 1e-4 for h in heat)

# 5
t = np.linspace(0, 2000, 201)
fo = analysis.fit_first_order(t, c0 * np.exp(-k * t))
cs = np.array([0.01, 0.03, 0.1, 0.3, 1.0])
o1 = analysis.reaction_order(cs, 0.5 * cs).slope
o2 = analysis.reaction_order(cs, 0.5 * cs**2).slope
Ts = np.array([288.15, 298.15, 308.15, 318.15])
Ea = 52.3
ks = 1e7 * np.exp(-Ea * 1000 / (analysis.R * Ts))
ar = analysis.arrhenius(Ts, ks).Ea_kJ
S = np.array([5, 10, 20, 50, 100, 200, 500.0])
mm = analysis.michaelis_menten(S, 3.2 * S / (41.0 + S))
fits = dict(k=fo.k / k - 1, order1=o1 - 1, order2=o2 / 2 - 1, Ea=ar / Ea - 1, Vmax=mm.Vmax / 3.2 - 1, Km=mm.Km / 41.0 - 1)
results["fits_ok"] = all(abs(v) <= 1e-6 for v in fits.values())

passed = all(results.values())
summary = dict(verdict="passed" if passed else "missed", **results, worst_exact_rel_err=worst_exact,
               worst_balance_rel_err=worst_balance, max_limit_rel_diff=max(abs(x["rel_diff"]) for x in conv),
               max_heat_rel_diff=max(abs(h["rel_diff"]) for h in heat),
               max_fit_rel_err=max(abs(v) for v in fits.values()), modes=modes, limit=conv, heat=heat, fits=fits)
print(checks.write_summary(HERE, summary))
print({k_: v for k_, v in summary.items() if not isinstance(v, (list, dict))})

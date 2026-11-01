"""A batch experiment: a stirred liquid at constant temperature with a gas space above it.

Three gas arrangements, as in the lab:

- ``burette``: the gas space connects to a gas burette held at the room's pressure; the collected volume grows as O2
  leaves the liquid. The gas is saturated with water vapor and sits at the vessel's temperature.
- ``closed``: a sealed vessel of fixed gas volume with a pressure gauge.
- ``open``: an open vessel; O2 leaving the liquid joins the room air, whose O2 partial pressure stays fixed.

The liquid starts saturated with air at the run's temperature and pressure. Dissolved O2 either stays in equilibrium
with the gas (``kLa_per_s=None``, the default) or crosses the surface at the rate kLa (c - c*), which lets the liquid
supersaturate as it does in a real vessel. Gases are treated as ideal (Benson and Krause's θ puts oxygen's departure at
about 0.07% at 1 atm and 25 °C), and the liquid as water for oxygen solubility.

The catalyst enters as a rate law; phase 1 offers a first-order rate constant that you supply.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Literal

import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp

from peroxide_engine import oxygen, sources, thermo, water

R = sources.value("gas_constant")
ATM = sources.value("standard_atmosphere")
T0 = sources.value("celsius_zero")
X_O2_AIR = sources.value("o2_fraction_dry_air")

# The range where every model the vessel uses has passed its check (see experiments/ and the Methodology page).
T_RANGE_K = (278.15, 328.15)
C0_MAX_M = 1.0
HEAT_T_K = thermo.T_REF


@dataclass
class Setup:
    c0_M: float  # initial H2O2, mol/L
    k_per_s: float  # first-order rate constant of the catalyst at the run temperature
    liquid_mL: float = 50.0
    gas_mL: float = 50.0  # gas space above the liquid at the start (headspace and tubing)
    T_K: float = 298.15
    pressure_Pa: float = ATM  # room pressure (burette, open) or starting pressure (closed)
    headspace: Literal["burette", "closed", "open"] = "burette"
    kLa_per_s: float | None = None  # None: dissolved O2 stays in equilibrium with the gas
    t_end_s: float | None = None  # default: five half-lives
    n_points: int = 401

    def validate(self) -> list[str]:
        out = []
        if not T_RANGE_K[0] <= self.T_K <= T_RANGE_K[1]:
            out.append(f"Temperature {self.T_K - T0:.1f} °C lies outside the checked range "
                       f"{T_RANGE_K[0] - T0:.0f}–{T_RANGE_K[1] - T0:.0f} °C.")
        if not 0 < self.c0_M <= C0_MAX_M:
            out.append(f"Starting H2O2 must lie in (0, {C0_MAX_M}] mol/L, the range where the solution behaves as "
                       "dilute for heat and oxygen solubility.")
        if self.k_per_s <= 0:
            out.append("The rate constant must be positive.")
        if self.liquid_mL <= 0 or self.gas_mL <= 0:
            out.append("Liquid and gas volumes must be positive.")
        if self.kLa_per_s is not None and self.kLa_per_s <= 0:
            out.append("kLa must be positive (or left empty for equilibrium).")
        return out


@dataclass
class Run:
    setup: Setup
    series: pd.DataFrame
    balance: pd.DataFrame | None = None
    readouts: list[dict] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)


def _gas(s: Setup):
    """Starting gas inventory: moles of O2 and of the other gases of air in the gas space, water vapor pressure."""
    T, P = s.T_K, s.pressure_Pa
    pw = float(water.vapor_pressure(T))
    n_dry = (P - pw) * s.gas_mL * 1e-6 / (R * T)
    return X_O2_AIR * n_dry, (1 - X_O2_AIR) * n_dry, pw


def simulate(s: Setup, rate: Callable[[np.ndarray], np.ndarray] | None = None) -> Run:
    """Run the experiment. ``rate`` maps H2O2 concentration (mol/L) to its rate of loss (mol/L/s); by default the
    first-order law k c with ``s.k_per_s``."""
    problems = s.validate()
    if problems:
        raise ValueError(" ".join(problems))
    T, VL = s.T_K, s.liquid_mL / 1000.0  # liquid in L
    custom_rate = rate is not None
    rate = rate or (lambda c: s.k_per_s * c)
    t_end = s.t_end_s or 5 * np.log(2) / s.k_per_s
    t = np.linspace(0.0, t_end, s.n_points)

    nO2_g0, n_inert, pw = _gas(s)
    P0 = s.pressure_Pa
    p_air_O2 = X_O2_AIR * (P0 - pw)
    cL0 = float(oxygen.concentration(T, p_air_O2, P0))
    nO2_total0 = nO2_g0 + cL0 * VL
    VG_m3 = s.gas_mL * 1e-6

    def henry(P):  # mol L-1 Pa-1
        return float(oxygen.henry(T, P))

    def equilibrium(n_tot):
        """Split O2 between liquid and gas for a total n_tot (mol); returns (nG, cL, P_total, dnG/dn_tot)."""
        if s.headspace == "open":
            return n_tot - cL0 * VL, cL0, P0, 1.0
        if s.headspace == "burette":
            a = henry(P0) * (P0 - pw) * VL
            b = n_inert + a - n_tot
            disc = np.sqrt(b * b + 4 * n_tot * n_inert)
            nG = 2 * n_tot * n_inert / (b + disc) if b > 0 else (disc - b) / 2
            cL = (n_tot - nG) / VL
            return nG, cL, P0, (nG + n_inert) / (2 * nG + n_inert + a - n_tot)
        P = P0
        for _ in range(6):  # closed: solubility depends weakly on total pressure through θ
            f = 1 + henry(P) * R * T * VL / VG_m3
            nG = n_tot / f
            P = (n_inert + nG) * R * T / VG_m3 + pw
        return nG, (n_tot - nG) / VL, P, 1 / f

    def p_O2(nG):
        if s.headspace == "open":
            return p_air_O2, P0
        if s.headspace == "burette":
            return (P0 - pw) * nG / (n_inert + nG), P0
        P = (n_inert + nG) * R * T / VG_m3 + pw
        return nG * R * T / VG_m3, P

    if s.kLa_per_s is None:
        if custom_rate:
            sol = solve_ivp(lambda _, y: [-rate(y[0])], (0, t_end), [s.c0_M], t_eval=t, method="Radau",
                            rtol=1e-10, atol=s.c0_M * 1e-13)
            if not sol.success:
                raise RuntimeError(sol.message)
            c = sol.y[0]
        else:
            c = s.c0_M * np.exp(-s.k_per_s * t)  # exact solution of the first-order law
        n_tot = nO2_total0 + 0.5 * VL * (s.c0_M - c)
        split = [equilibrium(n) for n in n_tot]
        nG = np.array([x[0] for x in split])
        cL = np.array([x[1] for x in split])
        P = np.array([x[2] for x in split])
        r = rate(c)
        J = np.array([x[3] for x in split]) * 0.5 * r  # O2 leaving the liquid, mol/L/s
    else:
        kla = s.kLa_per_s

        def rhs(_, y):
            c, cl, ng = y
            pO2, Ptot = p_O2(ng)
            j = kla * (cl - henry(Ptot) * pO2)
            rc = rate(c)
            return [-rc, 0.5 * rc - j, j * VL]

        y0 = [s.c0_M, cL0, nO2_g0 if s.headspace != "open" else 0.0]
        scale = max(s.c0_M, cL0)
        sol = solve_ivp(rhs, (0, t_end), y0, t_eval=t, method="Radau", rtol=1e-10,
                        atol=[scale * 1e-13, scale * 1e-13, scale * VL * 1e-13])
        if not sol.success:
            raise RuntimeError(sol.message)
        c, cL, nG = sol.y
        pairs = [p_O2(n) for n in nG]
        P = np.array([p[1] for p in pairs])
        r = rate(c)
        J = kla * (cL - np.array([henry(p[1]) * p[0] for p in pairs]))

    evolved = 0.5 * VL * (s.c0_M - c)  # mol O2 formed
    data = {
        "time_s": t,
        "H2O2_mM": c * 1e3,
        "rate_mM_per_s": r * 1e3,
        "O2_evolved_mmol": evolved * 1e3,
        "dissolved_O2_mg_per_L": cL * oxygen.molar_mass() * 1e3,
        "dissolved_O2_pct_air_saturation": cL / cL0 * 100,
    }
    if s.headspace == "burette":
        VG = (n_inert + nG) * R * T / (P0 - pw) * 1e6
        data["gas_collected_mL"] = VG - VG[0]
    elif s.headspace == "closed":
        data["pressure_kPa"] = P / 1e3
    else:  # open: nG counts O2 released to the room when kLa is set; at equilibrium every formed O2 leaves
        data["O2_released_mmol"] = (nG if s.kLa_per_s is not None else evolved) * 1e3
    heat_ok = abs(T - HEAT_T_K) < 0.005 and s.headspace in ("burette", "closed")
    if heat_ok:
        dH_aq = -thermo.decomposition_enthalpy_to_dissolved_o2(T) * 1e3  # J per mol H2O2, O2 left dissolved
        dH_des = -float(oxygen.solution_enthalpy(T)) * 1e3  # J per mol O2 leaving the liquid (absorbed)
        dH_vap = (thermo.nbs("H2O(g)", "dfH") - thermo.nbs("H2O(l)", "dfH")) * 1e3
        carried = pw / (P0 - pw) if s.headspace == "burette" else 0.0  # water vapor that saturates the new gas
        q = VL * (r * dH_aq - J * (dH_des + carried * dH_vap))
        data["heat_flow_mW"] = q * 1e3
        data["heat_J"] = np.concatenate([[0.0], np.cumsum(0.5 * (q[1:] + q[:-1]) * np.diff(t))])
    run = Run(s, pd.DataFrame(data))
    room = s.headspace == "open"
    run.balance = pd.DataFrame({  # oxygen inventory, mol O2; start + formed = liquid + gas + released
        "time_s": t,
        "start": np.full_like(t, nO2_total0),
        "formed": evolved,
        "liquid": cL * VL,
        "gas": np.full_like(t, nO2_g0) if room else nG,
        "released": (nG if s.kLa_per_s is not None else evolved) if room else np.zeros_like(t),
    })
    run.readouts = _readouts(run, cL0, heat_ok)
    if s.headspace == "closed" and P.max() > 3 * ATM:
        run.notes.append("The pressure passes 3 atm; Benson and Krause fitted the oxygen non-ideality term near 1 atm.")
    return run


def _readouts(run: Run, cL0: float, heat_ok: bool) -> list[dict]:
    from peroxide_engine import analysis

    s, d = run.setup, run.series
    t, c = d["time_s"].to_numpy(), d["H2O2_mM"].to_numpy()
    fit = analysis.fit_first_order(t, c)
    init = analysis.initial_rate(t, c)
    out = [
        dict(name="Rate constant k_obs", value=fit.k, units="1/s",
             how="Slope of ln[H2O2] against time over the whole run, as a lab fits its data."),
        dict(name="Half-life", value=analysis.half_life(t, c), units="s",
             how="Time at which [H2O2] reaches half its starting value, read from the curve."),
        dict(name="Initial rate", value=init, units="mM/s",
             how="Straight-line slope over the first 10% of the reaction, as a lab measures it."),
        dict(name="O2 formed by the end", value=d["O2_evolved_mmol"].iloc[-1], units="mmol",
             how="Half a mole of O2 per mole of H2O2 decomposed."),
        dict(name="Highest dissolved O2", value=d["dissolved_O2_mg_per_L"].max(), units="mg/L",
             how="Benson-Krause-Peterson solubility in water; equals the gas-side value at equilibrium."),
        dict(name="Highest dissolved O2, share of air saturation", value=d["dissolved_O2_pct_air_saturation"].max(),
             units="%", how="Relative to water saturated with air at the run's temperature and pressure."),
    ]
    if "gas_collected_mL" in d:
        T, P0 = s.T_K, s.pressure_Pa
        pw = float(water.vapor_pressure(T))
        vol = d["gas_collected_mL"].iloc[-1]
        out.append(dict(name="Gas collected in the burette", value=vol, units="mL",
                        how=f"Wet gas at {T - T0:.1f} °C and {P0 / 1e3:.2f} kPa, water vapor {pw / 1e3:.3f} kPa "
                            "(IAPWS)."))
        out.append(dict(name="Gas collected, dry at 0 °C and 1 atm", value=vol * (P0 - pw) / ATM * T0 / T,
                        units="mL", how="The burette volume with water vapor removed, by the ideal gas law."))
    if "pressure_kPa" in d:
        p = d["pressure_kPa"]
        out.append(dict(name="Final pressure", value=p.iloc[-1], units="kPa", how="Ideal gas in the fixed gas space, "
                        "plus water vapor (IAPWS)."))
        out.append(dict(name="Pressure rise", value=p.iloc[-1] - p.iloc[0], units="kPa", how="Final minus start."))
    if heat_ok:
        out.append(dict(name="Heat released", value=d["heat_J"].iloc[-1], units="J",
                        how="NBS enthalpy of H2O2(aq) to water and dissolved O2, minus the heat taken as O2 leaves "
                            "the liquid and, at the burette, as water evaporates into the new gas."))
        out.append(dict(name="Peak heat flow", value=d["heat_flow_mW"].max(), units="mW",
                        how="What an isothermal calorimeter records at its peak."))
    return out

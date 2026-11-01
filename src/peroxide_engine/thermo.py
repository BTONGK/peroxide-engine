"""Thermochemistry of H2O2 -> H2O + 1/2 O2 from the NBS tables (Wagman et al. 1982).

All values hold at 298.15 K and 0.1 MPa. The heat released by a solution depends on its strength, because mixing
H2O2 with water releases heat; the NBS dilution series (H2O2 in n H2O) gives that directly. The enthalpy of formation
of the solution "H2O2 in n H2O" is formed from the elements and n moles of liquid water (Wagman et al. 1982,
p. 2-22), so decomposing it gives

    H2O2(in n H2O) -> (n + 1) H2O(l) + 1/2 O2(g),    ΔrH = ΔfH(H2O, l) - ΔfH(H2O2 in n H2O).

Uncertainty: each table value is uncertain by 8 to 80 units of its last digit (two-standard-deviation scale, §4.2).
The enthalpy of H2O2 was measured against water's, so a reaction between them carries the H2O2 row's inaccuracy
alone, as §4.2 advises; the engine reports the upper end, 80 units (0.8 kJ/mol for a two-decimal row).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.interpolate import PchipInterpolator

from peroxide_engine import oxygen, sources, water

T_REF = 298.15
_nbs = sources.table("nbs")
R = sources.value("gas_constant")


@dataclass(frozen=True)
class Quantity:
    value: float
    uncertainty: float  # upper end of the source's stated range, two-standard-deviation scale
    units: str
    note: str


def nbs(key: str, column: str) -> float:
    return _nbs.value(key, column)


def nbs_uncertainty(key: str, column: str) -> float:
    """Upper end of the NBS inaccuracy for a tabulated value: 80 units of its last printed digit."""
    return 80 * _nbs.last_digit(key, column)


def molar_mass_h2o2() -> float:
    return nbs("H2O2(l)", "M")


def _reaction(peroxide: str, o2: str = "O2(g)") -> dict[str, float]:
    products = {"H2O(l)": 1.0, o2: 0.5}
    out = {}
    for col in ("dfH", "dfG"):
        out[col] = sum(nu * nbs(k, col) for k, nu in products.items()) - nbs(peroxide, col)
    out["S"] = sum(nu * nbs(k, "S") for k, nu in products.items()) - nbs(peroxide, "S")
    return out


def decomposition(state: str = "aq") -> dict[str, Quantity]:
    """ΔrH, ΔrG (kJ/mol) and ΔrS (J/mol/K) at 298.15 K for H2O2 -> H2O(l) + 1/2 O2(g).

    state "aq": H2O2 at infinite dilution in water (the NBS "ao" standard state); "l": pure liquid H2O2.
    """
    key = {"aq": "H2O2(ao)", "l": "H2O2(l)"}[state]
    r = _reaction(key)
    note = f"NBS tables, {key} -> H2O(l) + 1/2 O2(g), 298.15 K, 0.1 MPa"
    return {
        "dH": Quantity(r["dfH"], nbs_uncertainty(key, "dfH"), "kJ/mol", note),
        "dG": Quantity(r["dfG"], nbs_uncertainty(key, "dfG"), "kJ/mol", note),
        "dS": Quantity(r["S"], nbs_uncertainty(key, "S"), "J/(mol K)", note),
    }


def equilibrium_log10K(state: str = "aq") -> Quantity:
    """log10 of the equilibrium constant at 298.15 K, from ΔrG = -RT ln K."""
    g = decomposition(state)["dG"]
    k = -g.value * 1000 / (R * T_REF * np.log(10))
    dk = g.uncertainty * 1000 / (R * T_REF * np.log(10))
    return Quantity(k, dk, "", "from ΔrG of the NBS tables")


def _dilution_points():
    """(u, ΔfH) for the NBS series, with u = n/(n+1) the water fraction of the solution (pure H2O2 at 0, infinite
    dilution at 1), so the whole series lies on a finite interval."""
    pts = [(0.0, nbs("H2O2(l)", "dfH"))]
    for row in _nbs.rows:
        n = getattr(row, "n_water", None)
        if n is not None:
            n = float(n)
            pts.append((n / (n + 1), _nbs.value(row.key, "dfH")))
    pts.append((1.0, nbs("H2O2(ao)", "dfH")))
    u, h = map(np.array, zip(*sorted(pts)))
    return u, h


_U, _H = _dilution_points()
_interp = PchipInterpolator(_U, _H)


def formation_enthalpy_in_water(n_water):
    """ΔfH (kJ/mol) of H2O2 dissolved in n_water moles of water per mole, interpolated in the NBS dilution series.

    The interpolation is shape-preserving (PCHIP) in u = n/(n+1); experiments/decomposition_heat/ measures its error by
    leaving each tabulated point out in turn.
    """
    n = np.asarray(n_water, dtype=float)
    return _interp(n / (n + 1))


def water_per_peroxide(mass_fraction):
    """Moles of water per mole of H2O2 in a solution with this mass fraction of H2O2."""
    w = np.asarray(mass_fraction, dtype=float)
    return (1 - w) / water.molar_mass() / (w / molar_mass_h2o2())


def decomposition_enthalpy_in_water(n_water):
    """ΔrH (kJ/mol of H2O2) for H2O2(in n H2O) -> (n+1) H2O(l) + 1/2 O2(g) at 298.15 K."""
    return nbs("H2O(l)", "dfH") - formation_enthalpy_in_water(n_water)


def decomposition_heat_of_solution(mass_fraction) -> Quantity:
    """Heat released (kJ per kg of solution) when a solution of this H2O2 mass fraction decomposes completely to
    liquid water and O2 gas at 25 °C. The path runs from the initial solution to the final water, so the heat of
    dilution as it weakens is included."""
    w = float(mass_fraction)
    n = water_per_peroxide(w) if w < 1 else 0.0
    dh = float(decomposition_enthalpy_in_water(n)) if w < 1 else decomposition("l")["dH"].value
    mol_per_kg = w * 1000 / molar_mass_h2o2()
    return Quantity(-dh * mol_per_kg, 0.8 * mol_per_kg, "kJ/kg",
                    "NBS dilution series; uncertainty is 0.8 kJ per mol of H2O2 (80 units of 0.01 kJ/mol)")


def decomposition_enthalpy_to_dissolved_o2(T=T_REF) -> float:
    """ΔrH (kJ/mol) for H2O2(aq) -> H2O(l) + 1/2 O2(aq): the gas-phase value plus half the enthalpy of solution of O2
    (from the temperature slope of the Benson-Krause-Peterson solubility)."""
    return decomposition("aq")["dH"].value + 0.5 * float(oxygen.solution_enthalpy(T))

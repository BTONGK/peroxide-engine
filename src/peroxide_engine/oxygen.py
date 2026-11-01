"""Oxygen in water: solubility at any partial pressure, and its enthalpy of solution.

Solubility follows Benson, Krause and Peterson (1979): the mole fraction of O2 at 101.325 kPa O2, scaled to the actual
partial pressure (Henry's law) with Benson and Krause's gas non-ideality factor (1 - θP)/(1 - θ), the same factor the
USGS memorandum applies to air. The liquid is treated as pure water.
"""

from __future__ import annotations

import numpy as np

from peroxide_engine import sources, water

ATM = sources.value("standard_atmosphere")
T0 = sources.value("celsius_zero")
R = sources.value("gas_constant")


def molar_mass() -> float:
    """Molar mass of O2 in g/mol (NBS tables)."""
    return sources.table("nbs").value("O2(g)", "M")


def mole_fraction_1atm(T):
    """Mole fraction of O2 in water at 101.325 kPa O2 partial pressure (Benson, Krause & Peterson 1979)."""
    rec = sources.get("o2_mole_fraction_bkp1979")
    T = np.asarray(T, dtype=float)
    if np.any(T < rec.valid["T_min_K"] - 1e-9) or np.any(T > rec.valid["T_max_K"] + 1e-9):
        raise ValueError(f"oxygen solubility holds for {rec.valid['T_min_K']}-{rec.valid['T_max_K']} K; got {T}")
    c = rec.coefficients
    return np.exp(-(c["A"] + c["B"] / T + c["C"] / T**2))


def theta(T):
    """Benson and Krause's θ for O2 (atm-1): the gas non-ideality term, minus B/RT."""
    c = sources.coefficients("o2_virial_theta")
    t = np.asarray(T, dtype=float) - T0
    return c["a"] + c["b"] * t + c["c"] * t**2


def concentration(T, p_O2_Pa, p_total_Pa=ATM):
    """Dissolved O2 (mol/L) in equilibrium with an O2 partial pressure p_O2_Pa under total pressure p_total_Pa."""
    th = theta(T)
    x = mole_fraction_1atm(T) * (np.asarray(p_O2_Pa, dtype=float) / ATM) * (1 - th * p_total_Pa / ATM) / (1 - th)
    molality = x / ((1 - x) * water.molar_mass() / 1000.0)  # mol per kg water
    return molality * water.liquid_density(T) / 1000.0  # mol per L of water


def henry(T, p_total_Pa=ATM):
    """Dissolved O2 per unit O2 partial pressure, mol L-1 Pa-1 (dilute, so linear in the partial pressure)."""
    return concentration(T, 1.0, p_total_Pa)


def air_saturation(T, p_total_Pa=ATM):
    """Dissolved O2 (mol/L) in water at equilibrium with moist air at total pressure p_total_Pa."""
    p_O2 = sources.value("o2_fraction_dry_air") * (p_total_Pa - water.vapor_pressure(T))
    return concentration(T, p_O2, p_total_Pa)


def solution_enthalpy(T):
    """Enthalpy of solution of O2 gas into water, kJ/mol, from the temperature slope of the solubility (van 't Hoff).

    -ln x = A + B/T + C/T^2 gives d ln x / d(1/T) = -(B + 2C/T), so ΔH = R (B + 2C/T).
    """
    c = sources.coefficients("o2_mole_fraction_bkp1979")
    T = np.asarray(T, dtype=float)
    return R * (c["B"] + 2 * c["C"] / T) / 1000.0

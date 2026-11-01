"""Water: vapor pressure and liquid density from the IAPWS saturation equations."""

from __future__ import annotations

import numpy as np

from peroxide_engine import sources


def _check_range(T, rec_id: str):
    valid = sources.get(rec_id).valid
    T = np.asarray(T, dtype=float)
    if np.any(T < valid["T_min_K"] - 1e-9) or np.any(T > valid["T_max_K"] + 1e-9):
        raise ValueError(f"{rec_id} holds for {valid['T_min_K']}-{valid['T_max_K']} K; got {T}")
    return T


def vapor_pressure(T):
    """Vapor pressure of water in Pa at temperature T in K (IAPWS SR1-86(1992), eq. 1)."""
    T = _check_range(T, "iapws_vapor_pressure")
    c = sources.coefficients("iapws_vapor_pressure")
    tau = 1.0 - T / c["Tc"]
    s = (c["a1"] * tau + c["a2"] * tau**1.5 + c["a3"] * tau**3 + c["a4"] * tau**3.5 + c["a5"] * tau**4
         + c["a6"] * tau**7.5)
    return c["pc"] * np.exp(c["Tc"] / T * s)


def liquid_density(T):
    """Density of liquid water in kg/m3 at temperature T in K (IAPWS SR1-86(1992), eq. 2)."""
    T = _check_range(T, "iapws_liquid_density")
    c = sources.coefficients("iapws_liquid_density")
    tau = 1.0 - T / c["Tc"]
    s = (1 + c["b1"] * tau ** (1 / 3) + c["b2"] * tau ** (2 / 3) + c["b3"] * tau ** (5 / 3)
         + c["b4"] * tau ** (16 / 3) + c["b5"] * tau ** (43 / 3) + c["b6"] * tau ** (110 / 3))
    return c["rhoc"] * s


def molar_mass() -> float:
    """Molar mass of water in g/mol (NBS tables)."""
    return sources.table("nbs").value("H2O(l)", "M")

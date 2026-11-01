"""The fits an experimentalist applies to lab data, applied the same way to computed curves.

Each returns the fitted numbers with their standard errors, so a computed experiment reports what a lab would extract,
including the biases of the method (an initial rate taken over the first 10% of a first-order decay falls about 5%
below k·c0, for example).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import curve_fit

from peroxide_engine import sources

R = sources.value("gas_constant")


@dataclass(frozen=True)
class Line:
    slope: float
    slope_se: float
    intercept: float
    intercept_se: float
    r2: float
    n: int


def line(x, y) -> Line:
    """Ordinary least-squares straight line with standard errors."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    n = len(x)
    if n < 2:
        raise ValueError("a line needs at least two points")
    X = np.column_stack([x, np.ones(n)])
    coef, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ coef
    dof = n - 2
    s2 = float(resid @ resid / dof) if dof > 0 else 0.0
    cov = s2 * np.linalg.inv(X.T @ X)
    ss_tot = float(((y - y.mean()) ** 2).sum())
    r2 = 1 - float(resid @ resid) / ss_tot if ss_tot > 0 else 1.0
    return Line(float(coef[0]), float(np.sqrt(cov[0, 0])), float(coef[1]), float(np.sqrt(cov[1, 1])), r2, n)


@dataclass(frozen=True)
class FirstOrder:
    k: float
    k_se: float
    c0: float
    r2: float
    n: int


def fit_first_order(t, c) -> FirstOrder:
    """k from the slope of ln c against t (points with c > 0)."""
    t, c = np.asarray(t, float), np.asarray(c, float)
    keep = c > 0
    fit = line(t[keep], np.log(c[keep]))
    return FirstOrder(-fit.slope, fit.slope_se, float(np.exp(fit.intercept)), fit.r2, fit.n)


def initial_rate(t, c, fraction: float = 0.10) -> float:
    """Initial rate (units of c per unit t) from a straight line through the points up to ``fraction`` conversion."""
    t, c = np.asarray(t, float), np.asarray(c, float)
    keep = c >= c[0] * (1 - fraction)
    if keep.sum() < 3:
        keep[:3] = True
    return -line(t[keep], c[keep]).slope


def half_life(t, c) -> float:
    """Time at which c first reaches half its starting value, by linear interpolation; nan if it stays above."""
    t, c = np.asarray(t, float), np.asarray(c, float)
    half = c[0] / 2
    below = np.nonzero(c <= half)[0]
    if len(below) == 0:
        return float("nan")
    i = below[0]
    if i == 0:
        return 0.0
    return float(t[i - 1] + (c[i - 1] - half) * (t[i] - t[i - 1]) / (c[i - 1] - c[i]))


def reaction_order(concentrations, rates) -> Line:
    """Order from the slope of log(rate) against log(concentration)."""
    return line(np.log(concentrations), np.log(rates))


@dataclass(frozen=True)
class Arrhenius:
    Ea_kJ: float
    Ea_se_kJ: float
    lnA: float
    r2: float
    n: int


def arrhenius(T_K, k) -> Arrhenius:
    """Apparent activation energy from ln k against 1/T."""
    fit = line(1.0 / np.asarray(T_K, float), np.log(k))
    return Arrhenius(-fit.slope * R / 1000, fit.slope_se * R / 1000, fit.intercept, fit.r2, fit.n)


@dataclass(frozen=True)
class MichaelisMenten:
    Vmax: float
    Vmax_se: float
    Km: float
    Km_se: float
    n: int


def michaelis_menten(S, v) -> MichaelisMenten:
    """Vmax and Km by nonlinear least squares on v = Vmax S / (Km + S)."""
    S, v = np.asarray(S, float), np.asarray(v, float)
    p0 = [v.max(), float(np.median(S))]
    popt, pcov = curve_fit(lambda s, vm, km: vm * s / (km + s), S, v, p0=p0, maxfev=20000)
    se = np.sqrt(np.diag(pcov))
    return MichaelisMenten(float(popt[0]), float(se[0]), float(popt[1]), float(se[1]), len(S))


def turnover_frequency(rate_mol_per_s: float, sites_mol: float) -> float:
    """Molecules of H2O2 converted per site per second."""
    return rate_mol_per_s / sites_mol


def turnover_number(converted_mol: float, sites_mol: float) -> float:
    """Molecules of H2O2 converted per site over the run."""
    return converted_mol / sites_mol

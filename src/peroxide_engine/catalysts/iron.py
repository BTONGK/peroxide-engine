"""Dissolved iron(III) in perchloric acid: the mechanism of De Laat & Gallard (1999), 25 °C, I = 0.1 M.

Fe(III) binds H2O2 as two hydroperoxy complexes, FeIII(HO2)2+ (I1) and FeIII(OH)(HO2)+ (I2), which fall apart to
Fe2+ and HO2 (k6). Fe2+ then carries a chain with H2O2 through OH and HO2/O2-:

    Fe2+ + H2O2 -> Fe3+ + OH + OH-          (k7)      OH + H2O2 -> HO2 + H2O          (k9)
    Fe(III) + HO2/O2- -> Fe2+ + O2          (k11a,b)  Fe2+ + HO2/O2- -> FeIII(HO2)2+  (k10a,b)

with the radical recombinations of their Table 2. Hydrolysis of Fe(III) and the complexes stay at equilibrium (they
form within seconds), and [H+] stays fixed, as in the experiments. Every step is a mass-action rate; the peroxide in
the complexes counts toward the total H2O2 a titration measures, so oxygen and iron balance exactly.

Units: mol/L and seconds.
"""

from __future__ import annotations

from dataclasses import dataclass, fields, replace

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from peroxide_engine import sources

T_K = 298.15
IONIC_STRENGTH_M = 0.1
# The conditions of De Laat & Gallard's runs, where the same-laboratory check tests the model.
RANGE = {"pH": (1.0, 3.2), "fe_M": (50e-6, 1e-3), "c0_M": (0.2e-3, 1.0), "ratio": (1.0, 5000.0)}


def log_hydrolysis(id_: str, ionic_strength_M: float) -> float:
    """log10 of a hydrolysis constant at ionic strength I (Gallard et al. 1999, Table 1)."""
    c = sources.coefficients(id_)
    s = np.sqrt(ionic_strength_M)
    return c["a"] + c["b"] * s / (1 + c["c"] * s) + c["d"] * ionic_strength_M


def hydrolysis(id_: str, ionic_strength_M: float) -> float:
    return 10.0 ** log_hydrolysis(id_, ionic_strength_M)


@dataclass(frozen=True)
class Constants:
    K1: float  # Fe3+ + H2O = FeOH2+ + H+, M
    K2: float  # Fe3+ + 2 H2O = Fe(OH)2+ + 2 H+, M^2
    K22: float  # 2 Fe3+ + 2 H2O = Fe2(OH)2 4+ + 2 H+, M
    KI1: float  # Fe3+ + H2O2 = I1 + H+
    KI2: float  # FeOH2+ + H2O2 = I2 + H+
    k6: float
    k7: float
    k8: float
    k9: float
    k10a: float
    k10b: float
    k11a: float
    k11b: float
    k12a: float
    k12b: float
    k13a: float
    k13b: float
    k14a: float
    k14b: float
    k15: float

    @classmethod
    def published(cls, **changes: float) -> Constants:
        """De Laat & Gallard's Table 2 at I = 0.1 M, with k11a = 1e3 (see the constants file)."""
        rates = {f.name: sources.value(f"fe_{f.name}") for f in fields(cls) if f.name.startswith("k")}
        base = cls(K1=hydrolysis("fe_hydrolysis_K11", IONIC_STRENGTH_M),
                   K2=hydrolysis("fe_hydrolysis_K12", IONIC_STRENGTH_M), K22=sources.value("fe_dimer_K22"),
                   KI1=sources.value("fe_KI1"), KI2=sources.value("fe_KI2"), **rates)
        return replace(base, **changes)


def h_from_pH(pH: float) -> float:
    """[H+] from a measured pH at I = 0.1 M, with the activity coefficient both papers used."""
    return 10.0 ** -pH / sources.value("h_activity_0p1M")


@dataclass(frozen=True)
class Speciation:
    fe3: float  # Fe3+ (aquo)
    feoh: float
    feoh2: float
    dimer: float  # Fe2(OH)2 4+, per dimer
    I1: float
    I2: float
    peroxide_free: float

    @property
    def uncomplexed(self) -> float:
        """Fe(III) outside the peroxide complexes (De Laat & Gallard eq 12); this is the Fe(III) that meets HO2/O2-."""
        return self.fe3 + self.feoh + self.feoh2 + 2 * self.dimer

    @property
    def protons_released(self) -> float:
        return self.feoh + 2 * self.feoh2 + 2 * self.dimer + self.I1 + 2 * self.I2


def speciate(fe3_total: float, peroxide_total: float, h: float, k: Constants) -> Speciation:
    """Fe(III) species at equilibrium with H2O2 at fixed [H+]. The complexes hold a little peroxide, so the free
    H2O2 is found by iteration (three passes reach 1e-12 relative whenever iron is below the peroxide)."""
    p = peroxide_total
    for _ in range(4):
        a = 1 + k.K1 / h + k.K2 / h**2 + p * (k.KI1 / h + k.K1 * k.KI2 / h**2)
        b = 2 * k.K22 / h**2
        x = 2 * fe3_total / (a + np.sqrt(a * a + 4 * b * fe3_total))
        I1, I2 = k.KI1 * x * p / h, k.K1 * k.KI2 * x * p / h**2
        p = max(peroxide_total - I1 - I2, 0.0)
    return Speciation(x, k.K1 * x / h, k.K2 * x / h**2, k.K22 * x * x / h**2, I1, I2, p)


def h_from_acid(acid_M: float, fe3_total: float, peroxide_total: float, k: Constants) -> float:
    """Free [H+] when a known amount of strong acid is added: acid plus the protons iron releases as it hydrolyses
    and binds peroxide."""
    def balance(h):
        return h - acid_M - speciate(fe3_total, peroxide_total, h, k).protons_released

    return brentq(balance, acid_M, acid_M + 4 * fe3_total + 1e-12, xtol=1e-16, rtol=1e-13)


def initiation_constant(h: float, peroxide: float, k: Constants) -> float:
    """ki of De Laat & Gallard eq 20 (dimer left aside, as in their eq 15): rate of Fe2+ formation per
    [Fe(III)]total [H2O2]."""
    d = h * h + k.K1 * h + k.K2
    n = k.KI1 * h + k.K1 * k.KI2
    return k.k6 * n / d * d / (d + peroxide * n)


def apparent_ho2_constants(h: float, k: Constants) -> tuple[float, float]:
    """k10app and k11app: HO2/O2- with Fe(II) and with Fe(III), weighted by the HO2/O2- split at this [H+]."""
    ka = k.k12a / k.k12b
    f_o2m = ka / (ka + h)
    return (1 - f_o2m) * k.k10a + f_o2m * k.k10b, (1 - f_o2m) * k.k11a + f_o2m * k.k11b


# state: total peroxide (free plus complexed), Fe2+, OH, HO2, O2-, O2 formed
NAMES = ("H2O2", "Fe2", "OH", "HO2", "O2m", "O2")
STEPS = ("6", "7", "8", "9", "10a", "10b", "11a", "11b", "12a", "12b", "13a", "13b", "14a", "14b", "15")
# Change of each state variable per occurrence of each step (rows NAMES, columns STEPS). Step 6 and steps 10a/b
# move peroxide in and out of the complexes, which count toward total H2O2.
S = np.array([
    # 6   7   8   9  10a 10b 11a 11b 12a 12b 13a 13b 14a 14b 15
    [-1, -1, 0, -1, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 1],  # H2O2
    [1, -1, -1, 0, -1, -1, 1, 1, 0, 0, 0, 0, 0, 0, 0],  # Fe2+
    [0, 1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, -2],  # OH
    [1, 0, 0, 1, -1, 0, -1, 0, -1, 1, -2, -1, -1, 0, 0],  # HO2
    [0, 0, 0, 0, 0, -1, 0, -1, 1, -1, 0, -1, 0, -1, 0],  # O2-
    [0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0],  # O2
], dtype=float)


def _iron_terms(P: float, fe2: float, fe_total: float, h: float, k: Constants):
    """Complexed Fe(III), uncomplexed Fe(III) and free H2O2, with their derivatives in total H2O2 and in Fe2+."""
    def at(p_, f_):
        sp = speciate(fe_total - f_, p_, h, k)
        return np.array([sp.I1 + sp.I2, sp.uncomplexed, sp.peroxide_free])

    v = at(P, fe2)
    dp = max(P * 1e-6, 1e-15)
    df = max(fe_total * 1e-6, 1e-15)
    return v, (at(P + dp, fe2) - v) / dp, (at(P, fe2 + df) - v) / df


def _rates_and_jacobian(y: np.ndarray, fe_total: float, h: float, k: Constants, jac: bool = False):
    P, fe2, oh, ho2, o2m, _ = np.maximum(y, 0.0)
    (cx, fu, p), dP, dF = _iron_terms(P, fe2, fe_total, h, k) if jac else (
        _iron_terms_value(P, fe2, fe_total, h, k), None, None)
    r = np.array([k.k6 * cx, k.k7 * fe2 * p, k.k8 * fe2 * oh, k.k9 * oh * p, k.k10a * fe2 * ho2,
                  k.k10b * fe2 * o2m, k.k11a * fu * ho2, k.k11b * fu * o2m, k.k12a * ho2, k.k12b * o2m * h,
                  k.k13a * ho2 * ho2, k.k13b * ho2 * o2m, k.k14a * oh * ho2, k.k14b * oh * o2m, k.k15 * oh * oh])
    if not jac:
        return r, None
    d = np.zeros((15, 6))  # d(rate)/d(state)
    d[0, 0], d[0, 1] = k.k6 * dP[0], k.k6 * dF[0]
    d[1, 0], d[1, 1] = k.k7 * fe2 * dP[2], k.k7 * (p + fe2 * dF[2])
    d[2, 1], d[2, 2] = k.k8 * oh, k.k8 * fe2
    d[3, 0], d[3, 1], d[3, 2] = k.k9 * oh * dP[2], k.k9 * oh * dF[2], k.k9 * p
    d[4, 1], d[4, 3] = k.k10a * ho2, k.k10a * fe2
    d[5, 1], d[5, 4] = k.k10b * o2m, k.k10b * fe2
    d[6, 0], d[6, 1], d[6, 3] = k.k11a * dP[1] * ho2, k.k11a * dF[1] * ho2, k.k11a * fu
    d[7, 0], d[7, 1], d[7, 4] = k.k11b * dP[1] * o2m, k.k11b * dF[1] * o2m, k.k11b * fu
    d[8, 3] = k.k12a
    d[9, 4] = k.k12b * h
    d[10, 3] = 2 * k.k13a * ho2
    d[11, 3], d[11, 4] = k.k13b * o2m, k.k13b * ho2
    d[12, 2], d[12, 3] = k.k14a * ho2, k.k14a * oh
    d[13, 2], d[13, 4] = k.k14b * o2m, k.k14b * oh
    d[14, 2] = 2 * k.k15 * oh
    return r, d


def _iron_terms_value(P, fe2, fe_total, h, k):
    sp = speciate(fe_total - fe2, P, h, k)
    return sp.I1 + sp.I2, sp.uncomplexed, sp.peroxide_free


def _rhs(y, fe_total, h, k):
    return S @ _rates_and_jacobian(y, fe_total, h, k)[0]


def _jac(y, fe_total, h, k):
    return S @ _rates_and_jacobian(y, fe_total, h, k, jac=True)[1]


def oxidation_level(y: np.ndarray, fe_total: float) -> np.ndarray:
    """Oxidising equivalents above water and Fe(II): 2 per H2O2, 1 per OH, 3 per HO2/O2-, 4 per O2, 1 per Fe(III).
    Every step of the mechanism conserves it."""
    P, fe2, oh, ho2, o2m, o2 = y
    return 2 * P + oh + 3 * (ho2 + o2m) + 4 * o2 + (fe_total - fe2)


def _scales(c0: float, fe_total: float) -> np.ndarray:
    """Typical size of each state variable, so the solver's tolerances act on comparable numbers."""
    return np.array([c0, max(fe_total, 1e-9) * 1e-3, 1e-12, 1e-8, 1e-10, c0])


@dataclass
class Result:
    t: np.ndarray
    y: np.ndarray  # rows as NAMES, mol/L
    fe_total: float
    h: float
    constants: Constants
    sol: object  # the solver's dense output, in scaled variables
    scale: np.ndarray

    def __getitem__(self, name: str) -> np.ndarray:
        return self.y[NAMES.index(name)]

    def at(self, t) -> np.ndarray:
        """State at any time within the run, mol/L, rows as NAMES."""
        return self.sol.sol(np.atleast_1d(t)) * self.scale[:, None]

    def rate(self, t=None) -> np.ndarray:
        """Rate of H2O2 loss, mol/L/s, at the given times (default: the result's times)."""
        ys = self.y if t is None else self.at(t)
        return np.array([-_rhs(ys[:, i], self.fe_total, self.h, self.constants)[0] for i in range(ys.shape[1])])


def simulate(c0_M: float, fe_total_M: float, h_M: float, t_end_s: float, k: Constants | None = None,
             fe2_0_M: float = 0.0, t_eval: np.ndarray | None = None, rtol: float = 1e-8,
             stop_at_removal: float | None = None) -> Result:
    """Solve the mechanism from the moment H2O2 meets the iron. ``fe_total_M`` is all iron, of which ``fe2_0_M``
    starts as Fe(II). With ``stop_at_removal`` the run ends once that share of the H2O2 is gone."""
    k = k or Constants.published()
    scale = _scales(c0_M, fe_total_M)
    y0 = np.array([c0_M, fe2_0_M, 0.0, 0.0, 0.0, 0.0]) / scale

    def f(_, z):
        return _rhs(z * scale, fe_total_M, h_M, k) / scale

    def jac(_, z):
        return _jac(z * scale, fe_total_M, h_M, k) * scale[None, :] / scale[:, None]

    events = None
    if stop_at_removal is not None:
        def removed(_, z):
            return z[0] - (1 - stop_at_removal)
        removed.terminal, removed.direction = True, -1
        events = [removed]
    sol = solve_ivp(f, (0.0, t_end_s), y0, method="Radau", jac=jac, rtol=rtol, atol=1e-10,
                    dense_output=True, events=events)
    if not sol.success:
        raise RuntimeError(sol.message)
    t_stop = sol.t[-1]
    tt = np.linspace(0.0, t_stop, 401) if t_eval is None else np.asarray(t_eval)[np.asarray(t_eval) <= t_stop]
    return Result(tt, sol.sol(tt) * scale[:, None], fe_total_M, h_M, k, sol, scale)


def kobs_initial(c0_M: float, fe_total_M: float, h_M: float, k: Constants | None = None, removal: float = 0.2,
                 fe2_0_M: float = 0.0, n: int = 400) -> float:
    """Pseudo-first-order constant over the first ``removal`` of the H2O2: least-squares slope of ln(C/C0) against
    time through the origin, as De Laat & Gallard fitted their data (eq 13)."""
    k = k or Constants.published()
    guess = initiation_constant(h_M, c0_M, k) * fe_total_M  # initiation alone; the chain is faster
    t_max = 50 * removal / max(guess, 1e-15)
    res = simulate(c0_M, fe_total_M, h_M, t_max, k, fe2_0_M=fe2_0_M, stop_at_removal=removal)
    t_stop = res.sol.t[-1]
    if res.at(t_stop)[0, 0] > c0_M * (1 - removal) * (1 + 1e-6):
        raise RuntimeError("the run ended before reaching the requested removal")
    tt = np.linspace(0.0, t_stop, n + 1)[1:]
    yy = np.log(res.at(tt)[0] / c0_M)
    return float(-(tt @ yy) / (tt @ tt))


def validate(c0_M: float, fe_M: float, pH: float) -> list[str]:
    """Problems with a requested run, measured against the conditions of De Laat & Gallard's runs."""
    out = []
    lo, hi = RANGE["pH"]
    if not lo <= pH <= hi:
        out.append(f"pH {pH:g} lies outside {lo:g}-{hi:g}, the range of the measured runs (above pH 3.2 iron "
                   "precipitates).")
    lo, hi = RANGE["fe_M"]
    if not lo * (1 - 1e-9) <= fe_M <= hi * (1 + 1e-9):
        out.append(f"Fe(III) {fe_M * 1e6:g} uM lies outside {lo * 1e6:g}-{hi * 1e6:g} uM, the measured range.")
    lo, hi = RANGE["c0_M"]
    if not lo * (1 - 1e-9) <= c0_M <= hi * (1 + 1e-9):
        out.append(f"H2O2 {c0_M * 1e3:g} mM lies outside {lo * 1e3:g}-{hi * 1e3:g} mM, the measured range.")
    lo, hi = RANGE["ratio"]
    if fe_M > 0 and not lo * (1 - 1e-9) <= c0_M / fe_M <= hi * (1 + 1e-9):
        out.append(f"The ratio H2O2/Fe(III) = {c0_M / fe_M:.3g} lies outside {lo:g}-{hi:g}, the measured range.")
    return out


@dataclass
class Outcome:
    """A run of the iron model for the app: the solved mechanism and the numbers a lab would report."""
    result: Result
    t_end_s: float
    kobs: float  # over the first 20% removal, De Laat & Gallard's definition
    kobs_band: tuple[float, float]  # with k11a at the ends of its published bound
    h: float
    complexed_share: float  # Fe(III) in hydroperoxy complexes at the start

    def course(self, t: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        y = self.result.at(np.asarray(t, dtype=float))
        return y[0], self.result.rate(np.asarray(t, dtype=float))


def run(c0_M: float, fe_M: float, pH: float, until_removal: float = 0.9) -> Outcome:
    """Solve the mechanism until ``until_removal`` of the H2O2 is gone, with kobs and its k11a band."""
    problems = validate(c0_M, fe_M, pH)
    if problems:
        raise ValueError(" ".join(problems))
    k = Constants.published()
    h = h_from_pH(pH)
    kobs = kobs_initial(c0_M, fe_M, h, k)
    band = tuple(sorted(kobs_initial(c0_M, fe_M, h, Constants.published(k11a=v)) for v in (0.0, 2.0e3)))
    res = simulate(c0_M, fe_M, h, 200 * until_removal / kobs, k, stop_at_removal=until_removal)
    t_end = float(res.sol.t[-1])
    sp = speciate(fe_M, c0_M, h, k)
    return Outcome(res, t_end, kobs, band, h, (sp.I1 + sp.I2) / fe_M)

import numpy as np
import pytest

from peroxide_engine import analysis


def test_first_order_fit_and_half_life():
    t = np.linspace(0, 3000, 301)
    c = 50 * np.exp(-1e-3 * t)
    fit = analysis.fit_first_order(t, c)
    assert fit.k == pytest.approx(1e-3, rel=1e-9)
    assert analysis.half_life(t, c) == pytest.approx(np.log(2) / 1e-3, rel=1e-4)


def test_initial_rate_over_the_first_tenth_reads_about_five_percent_low():
    t = np.linspace(0, 5 * np.log(2) / 1e-3, 401)
    r = analysis.initial_rate(t, 100 * np.exp(-1e-3 * t))
    assert 0.94 < r / 0.1 < 0.96


def test_order_arrhenius_and_michaelis_menten_recover_known_values():
    cs = np.array([0.01, 0.1, 1.0])
    assert analysis.reaction_order(cs, 3 * cs**1.5).slope == pytest.approx(1.5)
    T = np.array([290.0, 300.0, 310.0])
    assert analysis.arrhenius(T, 1e8 * np.exp(-60000 / (analysis.R * T))).Ea_kJ == pytest.approx(60.0)
    S = np.array([1, 3, 10, 30, 100, 300.0])
    mm = analysis.michaelis_menten(S, 7 * S / (25 + S))
    assert (mm.Vmax, mm.Km) == (pytest.approx(7.0), pytest.approx(25.0))
    assert analysis.turnover_frequency(2e-6, 1e-9) == pytest.approx(2000)

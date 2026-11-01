import numpy as np
import pytest

from peroxide_engine import vessel

MODES = [(hs, kla) for hs in ("burette", "closed", "open") for kla in (None, 0.05)]


@pytest.mark.parametrize("headspace,kla", MODES)
def test_oxygen_is_conserved_and_h2o2_follows_the_exact_solution(headspace, kla):
    s = vessel.Setup(c0_M=0.2, k_per_s=2e-3, headspace=headspace, kLa_per_s=kla)
    run = vessel.simulate(s)
    t = run.series["time_s"].to_numpy()
    assert np.allclose(run.series["H2O2_mM"], 200 * np.exp(-2e-3 * t), rtol=1e-7)
    b = run.balance
    total = b["start"] + b["formed"]
    assert np.max(np.abs(b["liquid"] + b["gas"] + b["released"] - total) / total) < 1e-9


def test_burette_volume_matches_ideal_gas_for_the_oxygen_that_left():
    s = vessel.Setup(c0_M=0.1, k_per_s=1e-3, headspace="burette", t_end_s=50000)
    run = vessel.simulate(s)
    b = run.balance
    moved = b["gas"].iloc[-1] - b["gas"].iloc[0]
    from peroxide_engine import water
    pw = float(water.vapor_pressure(298.15))
    expect = moved * vessel.R * 298.15 / (s.pressure_Pa - pw) * 1e6
    assert run.series["gas_collected_mL"].iloc[-1] == pytest.approx(expect, rel=1e-9)
    assert b["formed"].iloc[-1] == pytest.approx(0.0025, rel=1e-6)  # 50 mL x 0.1 M x 1/2


def test_fast_transfer_reaches_equilibrium():
    eq = vessel.simulate(vessel.Setup(c0_M=0.2, k_per_s=2e-3, headspace="closed"))
    fast = vessel.simulate(vessel.Setup(c0_M=0.2, k_per_s=2e-3, headspace="closed", kLa_per_s=2.0))
    assert fast.series["pressure_kPa"].iloc[-1] == pytest.approx(eq.series["pressure_kPa"].iloc[-1], rel=1e-3)


def test_heat_appears_only_where_its_accounting_is_complete():
    at25 = vessel.simulate(vessel.Setup(c0_M=0.1, k_per_s=1e-3, headspace="closed"))
    at30 = vessel.simulate(vessel.Setup(c0_M=0.1, k_per_s=1e-3, headspace="closed", T_K=303.15))
    opened = vessel.simulate(vessel.Setup(c0_M=0.1, k_per_s=1e-3, headspace="open"))
    assert "heat_J" in at25.series and "heat_J" not in at30.series and "heat_J" not in opened.series
    per_mol = at25.series["heat_J"].iloc[-1] / (at25.balance["formed"].iloc[-1] * 2) / 1000
    assert 94.0 < per_mol < 100.7


def test_setup_outside_the_checked_range_is_refused():
    for bad in (dict(T_K=340.0), dict(c0_M=2.0), dict(k_per_s=-1.0)):
        s = vessel.Setup(**{"c0_M": 0.1, "k_per_s": 1e-3, **bad})
        assert s.validate()
        with pytest.raises(ValueError):
            vessel.simulate(s)


def test_a_rate_law_passed_as_a_function():
    s = vessel.Setup(c0_M=0.1, k_per_s=1e-3)
    second = vessel.simulate(s, rate=lambda c: 0.01 * c**2)
    c = second.series["H2O2_mM"].to_numpy() / 1000
    t = second.series["time_s"].to_numpy()
    assert np.allclose(c, 0.1 / (1 + 0.01 * 0.1 * t), rtol=1e-7)

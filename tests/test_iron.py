"""The dissolved-iron mechanism: balances, equilibria, the solver, and the coupling to the vessel."""

import numpy as np
import pytest

from peroxide_engine import vessel
from peroxide_engine.catalysts import iron

K = iron.Constants.published()


def test_every_step_conserves_oxidising_equivalents_and_iron():
    fe = 300e-6
    res = iron.simulate(30e-3, fe, iron.h_from_pH(2.5), 2e4, K)
    L = iron.oxidation_level(res.y, fe)
    assert np.max(np.abs(L - L[0])) / L[0] < 1e-10
    assert np.all(res["Fe2"] >= -1e-15) and np.all(res["Fe2"] <= fe)
    # each column of the stoichiometry conserves the oxidation level by itself
    weights = np.array([2, -1, 1, 3, 3, 4])  # H2O2, Fe2+ (one less Fe(III)), OH, HO2, O2-, O2
    assert np.allclose(weights @ iron.S, 0)


def test_speciation_accounts_for_all_iron_and_peroxide():
    fe, p, h = 1e-3, 0.5, iron.h_from_pH(3.0)
    sp = iron.speciate(fe, p, h, K)
    assert sp.uncomplexed + sp.I1 + sp.I2 == pytest.approx(fe, rel=1e-12)
    assert sp.peroxide_free + sp.I1 + sp.I2 == pytest.approx(p, rel=1e-12)


def test_free_acid_satisfies_the_proton_balance():
    acid, fe, p = 4.9e-3, 0.98e-3, 1.4
    h = iron.h_from_acid(acid, fe, p, K)
    assert h == pytest.approx(acid + iron.speciate(fe, p, h, K).protons_released, rel=1e-10)
    assert h > acid


def test_kobs_is_independent_of_solver_tolerance():
    args = (10e-3, 200e-6, iron.h_from_pH(2.64))
    a = iron.kobs_initial(*args, K)
    res = iron.simulate(*args, 3e4, K, rtol=1e-10, stop_at_removal=0.2)
    tt = np.linspace(0, res.sol.t[-1], 401)[1:]
    yy = np.log(res.at(tt)[0] / args[0])
    assert -(tt @ yy) / (tt @ tt) == pytest.approx(a, rel=1e-5)


def test_long_chain_limit_matches_the_papers_steady_state_formula():
    # De Laat & Gallard: for 50 < [H2O2]0/[Fe(III)]0 < 500 at pH 3 the full model and eqs 21-22 agree within 10%
    h = iron.h_from_pH(3.0)
    k10app, k11app = iron.apparent_ho2_constants(h, K)
    formula = 2 * np.sqrt(iron.initiation_constant(h, 50e-3, K) * K.k7 * k11app / k10app) * 200e-6
    assert iron.kobs_initial(50e-3, 200e-6, h, K) == pytest.approx(formula, rel=0.10)


def test_the_range_of_the_measured_runs_is_enforced():
    assert iron.validate(0.05, 200e-6, 3.0) == []
    assert iron.validate(0.05, 200e-6, 3.5)
    assert iron.validate(0.05, 2e-3, 3.0)
    assert iron.validate(2.0, 200e-6, 3.0)
    assert iron.validate(10e-3, 1e-3, 3.0) == []  # ratio 10
    assert iron.validate(0.2e-3, 1e-3, 3.0)  # ratio 0.2


def test_the_vessel_follows_the_mechanisms_course():
    out = iron.run(0.05, 200e-6, 3.0, until_removal=0.5)
    s = vessel.Setup(c0_M=0.05, T_K=iron.T_K, headspace="burette", t_end_s=out.t_end_s)
    run = vessel.simulate(s, course=out.course)
    c = run.series["H2O2_mM"].to_numpy() / 1000
    assert c[-1] == pytest.approx(0.025, rel=1e-6)
    assert run.series["O2_evolved_mmol"].iloc[-1] == pytest.approx(0.5 * 0.05 * 0.025 * 1e3, rel=1e-6)
    b = run.balance
    assert np.max(np.abs(b.start + b.formed - b.liquid - b.gas - b.released)) < 1e-15
    assert out.kobs_band[0] <= out.kobs <= out.kobs_band[1]


def test_run_page_computes_an_iron_run():
    from streamlit.testing.v1 import AppTest

    from peroxide_engine import checks

    page = checks.ROOT / "src" / "peroxide_engine" / "ui" / "views" / "run.py"
    at = AppTest.from_file(str(page), default_timeout=120).run()
    catalyst = next(r for r in at.radio if r.label == "Catalyst")
    catalyst.set_value(catalyst.options[1]).run()
    assert not at.exception, at.exception
    assert not at.error
    names = at.dataframe[0].value["Readout"].tolist()
    assert names[0].startswith("k_obs") and any(n.startswith("Second-order") for n in names)

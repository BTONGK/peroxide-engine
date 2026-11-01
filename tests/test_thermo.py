import numpy as np
import pytest

from peroxide_engine import oxygen, sources, thermo, water


def test_reaction_values_from_the_tables():
    aq, liq = thermo.decomposition("aq"), thermo.decomposition("l")
    assert aq["dH"].value == pytest.approx(-94.66, abs=1e-9)
    assert liq["dH"].value == pytest.approx(-98.05, abs=1e-9)
    assert aq["dG"].value == pytest.approx(-103.099, abs=1e-9)
    assert aq["dH"].uncertainty == pytest.approx(0.8)
    assert thermo.equilibrium_log10K().value == pytest.approx(18.062, abs=1e-3)


def test_gibbs_helmholtz_consistency_within_the_tables_uncertainty():
    for state in ("aq", "l"):
        d = thermo.decomposition(state)
        assert abs(d["dH"].value - 298.15 * d["dS"].value / 1000 - d["dG"].value) < 0.8


def test_dilution_series_passes_through_every_tabulated_point():
    nbs = sources.table("nbs")
    for row in nbs.rows:
        n = getattr(row, "n_water", None)
        if n is not None:
            assert float(thermo.formation_enthalpy_in_water(float(n))) == pytest.approx(nbs.value(row.key, "dfH"))
    assert float(thermo.formation_enthalpy_in_water(1e9)) == pytest.approx(-191.17, abs=1e-6)
    assert float(thermo.formation_enthalpy_in_water(0.0)) == pytest.approx(-187.78)


def test_heat_per_mole_rises_with_strength():
    w = np.linspace(0.01, 0.99, 50)
    per_mol = [thermo.decomposition_heat_of_solution(x).value / (x * 1000 / thermo.molar_mass_h2o2()) for x in w]
    assert np.all(np.diff(per_mol) > 0)
    assert per_mol[0] == pytest.approx(94.66, abs=0.05)


def test_water_and_oxygen_reference_points():
    assert float(water.vapor_pressure(298.15)) == pytest.approx(3169.8, abs=0.1)
    assert float(oxygen.air_saturation(298.15)) * oxygen.molar_mass() * 1000 == pytest.approx(8.26, abs=0.01)
    assert float(oxygen.solution_enthalpy(298.15)) == pytest.approx(-12.01, abs=0.01)


def test_ranges_are_enforced():
    with pytest.raises(ValueError):
        water.vapor_pressure(270.0)
    with pytest.raises(ValueError):
        oxygen.mole_fraction_1atm(380.0)

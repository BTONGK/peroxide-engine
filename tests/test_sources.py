"""Every constant carries its source, and every printed check value is reproduced."""

import subprocess

import pytest

from peroxide_engine import checks, oxygen, sources, water
from peroxide_engine.catalysts import iron


def test_every_record_names_its_source():
    reg = sources.registry()
    assert len(reg) >= 10
    for rec in reg.values():
        assert rec.source.doi or rec.source.url
        assert rec.source.locator
        if rec.kind == "table":
            assert all(row.locator for row in rec.rows)


def test_a_record_without_a_link_is_refused():
    with pytest.raises(ValueError):
        sources.Source(cite="Someone 2020", reference="An unnamed paper with ten chars", locator="p. 1")


FUNCS = {"iapws_vapor_pressure": water.vapor_pressure, "iapws_liquid_density": water.liquid_density,
         "o2_mole_fraction_bkp1979": oxygen.mole_fraction_1atm,
         "fe_hydrolysis_K11": lambda I: iron.hydrolysis("fe_hydrolysis_K11", I),
         "fe_hydrolysis_K12": lambda I: iron.hydrolysis("fe_hydrolysis_K12", I)}


def test_printed_check_values_are_reproduced():
    with_checks = [r for r in sources.registry().values() if r.kind == "correlation" and r.check_values]
    assert {r.id for r in with_checks} == set(FUNCS)
    for rec in with_checks:
        for cv in rec.check_values:
            arg = next(v for key, v in cv.items() if key != "value")  # the one input: T_K or I_M
            assert abs(float(FUNCS[rec.id](arg)) / cv["value"] - 1) <= rec.check_rel_tolerance, (rec.id, cv)


def test_last_digit_reads_the_printed_value():
    nbs = sources.table("nbs")
    assert nbs.last_digit("H2O2(ao)", "dfH") == pytest.approx(0.01)
    assert nbs.last_digit("H2O(l)", "dfH") == pytest.approx(0.001)
    assert nbs.last_digit("O2(ao)", "dfH") == pytest.approx(0.1)


def test_papers_and_cache_stay_out_of_git():
    try:
        files = subprocess.run(["git", "ls-files"], cwd=checks.ROOT, capture_output=True, text=True, check=True)
    except (FileNotFoundError, subprocess.CalledProcessError):
        pytest.skip("outside a git checkout")
    tracked = files.stdout.split()
    assert not [f for f in tracked if f.startswith("cache/") or "/cache/" in f or f.lower().endswith(".pdf")]

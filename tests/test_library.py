"""The measured catalyst library: every entry sourced, every class present, ranges enforced, and the Run page."""

import numpy as np
import pytest

from peroxide_engine import checks
from peroxide_engine.catalysts import library


def test_every_class_is_listed_and_every_entry_is_sourced():
    entries = library.entries().values()
    assert {e.cls for e in entries} == set(library.CLASSES)
    for e in entries:
        assert e.source.doi or e.source.url
        assert e.status == "awaiting" or e.law is not None


def test_ready_entries_hold_their_printed_line():
    zro2 = library.entries()["zro2_lousada"].law
    av = 7.5 / 50e-6  # the thesis's standard 7.5 m² in 50 mL
    assert zro2.k_obs(av, 298.15) == pytest.approx(2.39e-10 * av + 3e-5)
    lo, hi = zro2.k_band(av, 298.15)
    assert lo < 6.15e-5 < hi  # the printed k1 sits inside the band of the printed uncertainties
    # the activation energy carries k_obs exactly by its Arrhenius factor
    ratio = zro2.k_obs(av, 323.15) / zro2.k_obs(av, 298.15)
    assert ratio == pytest.approx(np.exp(-33e3 / library.R * (1 / 323.15 - 1 / 298.15)))


def test_ranges_are_enforced():
    hfo2 = library.entries()["hfo2_lousada"].law
    assert hfo2.check(1e5, 0.5e-3, 300.0) == []
    assert len(hfo2.check(1e6, 1e-3, 360.0)) == 3
    assert hfo2.area_per_volume(15.0) == pytest.approx(150000)  # 0.75 g in 50 mL at 10.0 m²/g


def test_run_page_computes_a_library_run():
    from streamlit.testing.v1 import AppTest

    page = checks.ROOT / "src" / "peroxide_engine" / "ui" / "views" / "run.py"
    at = AppTest.from_file(str(page), default_timeout=60).run()
    catalyst = next(r for r in at.radio if r.label == "Catalyst")
    catalyst.set_value(catalyst.options[2]).run()
    assert not at.exception, at.exception
    assert not at.error
    names = at.dataframe[0].value["Readout"].tolist()
    assert names[0] == "k_obs from the measured line"

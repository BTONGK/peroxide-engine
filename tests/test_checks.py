"""Each model check's summary holds the verdict its report shows, and the app pages load."""

from pathlib import Path

import pytest
from streamlit.testing.v1 import AppTest

from peroxide_engine import checks

VERDICTS = {"transcription": "passed", "decomposition_heat": "passed", "oxygen_solubility": "passed",
            "exact_solutions": "passed", "iron_transcription": "passed", "iron_same_lab": "passed",
            "iron_other_lab": "missed", "library_consistency": "passed",
            "oxide_descriptor_check": "missed"}


def test_check_verdicts():
    found = {c["folder"]: c["summary"]["verdict"] for c in checks.all_checks()}
    assert found == VERDICTS
    for c in checks.all_checks():
        assert (checks.EXPERIMENTS / c["folder"] / "PASS_RULE.md").exists()


VIEWS = Path(checks.ROOT) / "src" / "peroxide_engine" / "ui" / "views"


@pytest.mark.parametrize("page", sorted(p.name for p in VIEWS.glob("*.py")))
def test_page_runs(page):
    at = AppTest.from_file(str(VIEWS / page), default_timeout=60).run()
    assert not at.exception, at.exception

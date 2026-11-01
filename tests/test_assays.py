import pytest

from peroxide_engine import assays


@pytest.mark.parametrize("equation", [assays.PERMANGANATE, *assays.IODOMETRIC])
def test_equations_balance(equation):
    atoms, charge = assays.balance(equation)
    assert not atoms and charge == 0


def test_titres():
    assert assays.permanganate_titre_mL(0.1, 1.0, 0.02) == pytest.approx(2.0)
    assert assays.thiosulfate_titre_mL(0.1, 1.0, 0.1) == pytest.approx(2.0)

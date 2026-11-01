"""What titrations read, from computed H2O2 concentrations. These use stoichiometry alone; the tests check that each
equation balances in atoms and charge.

Spectroscopic assays (absorbance at 240 nm, triiodide at 353 nm, titanium peroxo at 407 nm) join once their molar
absorption coefficients are taken from the primary papers and checked.

Charges are written in brackets, e.g. MnO4[-], Mn[2+], S2O3[2-].
"""

from __future__ import annotations

import re
from collections import Counter

PERMANGANATE = "2 MnO4[-] + 5 H2O2 + 6 H[+] -> 2 Mn[2+] + 5 O2 + 8 H2O"
IODOMETRIC = ("H2O2 + 2 I[-] + 2 H[+] -> I2 + 2 H2O", "I2 + 2 S2O3[2-] -> 2 I[-] + S4O6[2-]")


def permanganate_titre_mL(c_H2O2_M, aliquot_mL, c_KMnO4_M):
    """Volume of KMnO4 solution that reaches the end point for an acidified aliquot: 2 MnO4- per 5 H2O2."""
    return 2 / 5 * c_H2O2_M * aliquot_mL / c_KMnO4_M


def thiosulfate_titre_mL(c_H2O2_M, aliquot_mL, c_thiosulfate_M):
    """Volume of thiosulfate for an aliquot treated with excess acidified iodide: 2 S2O3^2- per H2O2."""
    return 2 * c_H2O2_M * aliquot_mL / c_thiosulfate_M


_SPECIES = re.compile(r"^(\d*)\s*([A-Za-z0-9()]+)(?:\[(\d*)([+-])\])?$")
_ELEMENT = re.compile(r"([A-Z][a-z]?)(\d*)")


def _side(text: str) -> tuple[Counter, int]:
    atoms, charge = Counter(), 0
    for term in text.split(" + "):
        m = _SPECIES.match(term.strip())
        if m is None:
            raise ValueError(f"unreadable species {term!r}")
        nu = int(m.group(1) or 1)
        for el, count in _ELEMENT.findall(m.group(2)):
            atoms[el] += nu * int(count or 1)
        if m.group(4):
            charge += nu * int(m.group(3) or 1) * (1 if m.group(4) == "+" else -1)
    return atoms, charge


def balance(equation: str) -> tuple[Counter, int]:
    """Atoms and charge of the products minus the reactants; zero for a balanced equation."""
    left, right = equation.split(" -> ")
    a1, q1 = _side(left)
    a2, q2 = _side(right)
    diff = Counter(a2)
    diff.subtract(a1)
    return Counter({k: v for k, v in diff.items() if v}), q2 - q1

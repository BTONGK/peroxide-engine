"""Shared helpers for the app pages."""

from __future__ import annotations

import math
from pathlib import Path

import streamlit as st
import yaml

ROOT = Path(__file__).resolve().parents[3]
DOCS = ROOT / "docs"


@st.cache_data
def doc(name: str) -> dict:
    return yaml.safe_load((DOCS / f"{name}.yaml").read_text())


def sig(x: float, digits: int = 4) -> str:
    """A number to `digits` significant figures, in plain or scientific notation as fits."""
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return "–"
    if x == 0:
        return "0"
    mag = math.floor(math.log10(abs(x)))
    if -3 <= mag < 6:
        return f"{x:.{max(digits - 1 - mag, 0)}f}"
    return f"{x:.{digits - 1}e}"


def verdict_badge(verdict: str) -> str:
    return {"passed": ":green-badge[:material/check: Passed]",
            "missed": ":red-badge[:material/close: Missed]"}.get(verdict, f":gray-badge[{verdict}]")

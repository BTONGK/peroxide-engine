"""The measured catalyst library: catalysts whose rate a laboratory measured, run inside the conditions it measured.

Entries live in YAML under ``data/library/``, one file per class. A ready entry carries its rate law with every value
as printed and the ranges of H2O2, temperature and catalyst amount the source measured; an awaiting entry names the
paper it needs. The app runs a ready entry only inside its ranges.
"""

from __future__ import annotations

from functools import cache
from pathlib import Path
from typing import Literal

import numpy as np
import yaml
from pydantic import BaseModel, ConfigDict, Field, model_validator

from peroxide_engine.sources import Source

DATA = Path(__file__).parent.parent / "data" / "library"
R = 8.314462618  # J/(mol K)
CLASSES = {"enzyme": "Enzymes", "dissolved": "Dissolved catalysts", "oxide": "Solid oxides", "metal": "Metals"}


class Value(BaseModel):
    """One number as its source prints it."""

    model_config = ConfigDict(extra="forbid")
    value: float
    uncertainty: float = Field(0.0, ge=0)
    units: str
    printed: str
    locator: str
    note: str | None = None


class SurfaceFirstOrder(BaseModel):
    """First order in H2O2 on a solid: k_obs = (k2 × Sa/V + b) × exp(−Ea/R (1/T − 1/T_ref)).

    k2 (m/s) and the intercept b (1/s) are the source's straight line of k_obs against the surface area per volume
    Sa/V (1/m) at T_ref; the activation energy carries it across the measured temperatures.
    """

    model_config = ConfigDict(extra="forbid")
    kind: Literal["surface_first_order"]
    k2: Value
    intercept: Value
    T_ref_K: float
    Ea: Value | None = None
    A: Value | None = None
    specific_area: Value | None = None
    T_range_K: tuple[float, float]
    area_per_volume_per_m: tuple[float, float]
    c0_M: tuple[float, float]
    medium: str
    ranges_from: str = Field(min_length=20)

    def k_obs(self, area_per_volume: float, T_K: float, k2: float | None = None, b: float | None = None,
              Ea_kJ: float | None = None) -> float:
        k2 = self.k2.value if k2 is None else k2
        b = self.intercept.value if b is None else b
        k = k2 * area_per_volume + b
        if self.Ea is not None and T_K != self.T_ref_K:
            Ea = (self.Ea.value if Ea_kJ is None else Ea_kJ) * 1e3
            k *= np.exp(-Ea / R * (1 / T_K - 1 / self.T_ref_K))
        return float(k)

    def k_band(self, area_per_volume: float, T_K: float) -> tuple[float, float]:
        """k_obs at the ends of the printed uncertainties of k2, b and Ea."""
        corners = [self.k_obs(area_per_volume, T_K, self.k2.value + s2 * self.k2.uncertainty,
                              self.intercept.value + sb * self.intercept.uncertainty,
                              None if self.Ea is None else self.Ea.value + se * self.Ea.uncertainty)
                   for s2 in (-1, 1) for sb in (-1, 1) for se in (-1, 1)]
        return min(corners), max(corners)

    def area_per_volume(self, mass_g_per_L: float) -> float:
        """Sa/V in 1/m from the catalyst loading, through the measured specific surface area."""
        return mass_g_per_L * self.specific_area.value * 1e3  # m²/g × g/L × 1000 L/m³

    def check(self, area_per_volume: float, c0_M: float, T_K: float) -> list[str]:
        out = []
        lo, hi = self.area_per_volume_per_m
        if not lo * (1 - 1e-9) <= area_per_volume <= hi * (1 + 1e-9):
            out.append(f"The surface area per volume is {area_per_volume:,.0f} per m; the source measured "
                       f"{lo:,.0f}–{hi:,.0f} per m.")
        lo, hi = self.c0_M
        if not lo * (1 - 1e-9) <= c0_M <= hi * (1 + 1e-9):
            out.append(f"Starting H2O2 is {c0_M * 1e3:g} mM; the source measured {lo * 1e3:g}–{hi * 1e3:g} mM.")
        lo, hi = self.T_range_K
        if not lo - 1e-6 <= T_K <= hi + 1e-6:
            out.append(f"The temperature is {T_K - 273.15:.2f} °C; the source measured {lo - 273.15:.2f}–"
                       f"{hi - 273.15:.2f} °C.")
        return out


class Mechanism(BaseModel):
    """A checked mechanism with its own module and Run page option."""

    model_config = ConfigDict(extra="forbid")
    kind: Literal["mechanism"]
    module: str
    summary: str


class Entry(BaseModel):
    model_config = ConfigDict(extra="forbid", populate_by_name=True)
    id: str
    name: str
    cls: Literal["enzyme", "dissolved", "oxide", "metal"] = Field(alias="class")
    material: str
    status: Literal["ready", "awaiting"]
    needs: str | None = None
    law: SurfaceFirstOrder | Mechanism | None = Field(None, discriminator="kind")
    source: Source
    notes: str | None = None

    @model_validator(mode="after")
    def _complete(self):
        if self.status == "ready" and self.law is None:
            raise ValueError(f"{self.id}: a ready entry needs its rate law")
        if self.status == "awaiting" and not self.needs:
            raise ValueError(f"{self.id}: an awaiting entry names what it needs")
        return self


@cache
def entries() -> dict[str, Entry]:
    out: dict[str, Entry] = {}
    for f in sorted(DATA.glob("*.yaml")):
        for raw in yaml.safe_load(f.read_text()):
            e = Entry.model_validate(raw)
            if e.id in out:
                raise ValueError(f"duplicate library id {e.id}")
            out[e.id] = e
    return out


def runnable() -> list[Entry]:
    """Ready entries the Run page computes through the library (mechanisms have their own option)."""
    return [e for e in entries().values() if e.status == "ready" and isinstance(e.law, SurfaceFirstOrder)]

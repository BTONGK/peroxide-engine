"""Every number the engine computes with, and where it comes from.

Records live in YAML under ``data/constants/``, one file per source family. Three kinds:

- ``constant``: one value with units, an uncertainty and how that uncertainty was obtained, the conditions it holds
  for, and its source down to the table, page or equation.
- ``table``: rows of values from one source. Values are stored as printed (strings), so the last printed digit, which
  carries the source's uncertainty statement, survives loading.
- ``correlation``: an equation from a source with its coefficients as printed, its valid range, the source's accuracy
  statement, and the check values the source prints so a test can confirm the transcription.

Every record names a source with a DOI or URL and a locator; the loader refuses a record without one.
"""

from __future__ import annotations

from decimal import Decimal
from functools import cache
from pathlib import Path
from typing import Literal

import yaml
from pydantic import BaseModel, ConfigDict, Field, model_validator

DATA = Path(__file__).parent / "data" / "constants"


class Source(BaseModel):
    model_config = ConfigDict(extra="forbid")
    cite: str = Field(min_length=3)
    reference: str = Field(min_length=10)
    doi: str | None = None
    url: str | None = None
    locator: str = Field(min_length=2)

    @model_validator(mode="after")
    def _link(self):
        if not (self.doi or self.url):
            raise ValueError(f"source {self.cite!r} needs a DOI or a URL")
        return self


class Conditions(BaseModel):
    model_config = ConfigDict(extra="forbid")
    T_K: float | None = None
    p_Pa: float | None = None
    pH: float | None = None
    ionic_strength_M: float | None = None
    medium: str | None = None


class Constant(BaseModel):
    model_config = ConfigDict(extra="forbid")
    kind: Literal["constant"]
    id: str
    name: str
    value: float
    units: str
    uncertainty: float = Field(ge=0)
    uncertainty_note: str = Field(min_length=10)
    exact: bool = False
    conditions: Conditions = Conditions()
    source: Source


class Correlation(BaseModel):
    model_config = ConfigDict(extra="forbid")
    kind: Literal["correlation"]
    id: str
    name: str
    equation: str
    coefficients: dict[str, float]
    units: str
    valid: dict[str, float]
    accuracy_note: str = Field(min_length=10)
    transcribed_from: str = Field(min_length=10, description="where the coefficients were read, e.g. a page image")
    check_values: list[dict[str, float]] = []
    check_rel_tolerance: float | None = Field(None, gt=0, description="largest |computed/printed - 1| allowed")
    conditions: Conditions = Conditions()
    source: Source

    @model_validator(mode="after")
    def _tolerance(self):
        if self.check_values and self.check_rel_tolerance is None:
            raise ValueError(f"{self.id}: check values need check_rel_tolerance")
        return self


class Row(BaseModel):
    model_config = ConfigDict(extra="allow")
    key: str
    locator: str


class Table(BaseModel):
    model_config = ConfigDict(extra="forbid")
    kind: Literal["table"]
    id: str
    name: str
    columns: dict[str, str]
    conditions: Conditions = Conditions()
    uncertainty_note: str = Field(min_length=10)
    source: Source
    rows: list[Row] = Field(min_length=1)

    def row(self, key: str) -> Row:
        for r in self.rows:
            if r.key == key:
                return r
        raise KeyError(f"{key} is missing from table {self.id}")

    def value(self, key: str, column: str) -> float:
        return float(Decimal(self.printed(key, column)))

    def printed(self, key: str, column: str) -> str:
        text = getattr(self.row(key), column, None)
        if text is None or text == "-":
            raise KeyError(f"table {self.id} lacks {column} for {key}")
        return str(text)

    def last_digit(self, key: str, column: str) -> float:
        """One unit of the last printed digit of a value, e.g. 0.01 for -191.17."""
        return float(Decimal(1).scaleb(Decimal(self.printed(key, column)).as_tuple().exponent))


Record = Constant | Correlation | Table


def _parse(item: dict, path: Path) -> Record:
    kind = item.get("kind")
    model = {"constant": Constant, "correlation": Correlation, "table": Table}.get(kind)
    if model is None:
        raise ValueError(f"{path.name}: record {item.get('id')!r} has kind {kind!r}")
    return model.model_validate(item)


@cache
def registry() -> dict[str, Record]:
    records: dict[str, Record] = {}
    for path in sorted(DATA.glob("*.yaml")):
        for item in yaml.safe_load(path.read_text()):
            rec = _parse(item, path)
            if rec.id in records:
                raise ValueError(f"{path.name}: duplicate id {rec.id}")
            records[rec.id] = rec
    return records


def get(id_: str) -> Record:
    return registry()[id_]


def value(id_: str) -> float:
    rec = get(id_)
    assert isinstance(rec, Constant), f"{id_} is a {rec.kind}"
    return rec.value


def coefficients(id_: str) -> dict[str, float]:
    rec = get(id_)
    assert isinstance(rec, Correlation), f"{id_} is a {rec.kind}"
    return rec.coefficients


def table(id_: str) -> Table:
    rec = get(id_)
    assert isinstance(rec, Table), f"{id_} is a {rec.kind}"
    return rec

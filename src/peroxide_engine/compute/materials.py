"""Crystal structures from the Materials Project (materialsproject.org).

The API key comes from MP_API_KEY, PMG_MAPI_KEY or ~/.config/.pmgrc.yaml (the pymatgen convention) and travels only in
the X-API-KEY header, to api.materialsproject.org. Structures are cached in the git-ignored cache/materials/.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import requests
import yaml

API = "https://api.materialsproject.org/materials/summary/"
FIELDS = "material_id,formula_pretty,energy_above_hull,symmetry,nsites,band_gap,theoretical"
CACHE = Path(__file__).resolve().parents[3] / "cache" / "materials"


class MaterialsProjectError(RuntimeError):
    pass


def api_key() -> str | None:
    for var in ("MP_API_KEY", "PMG_MAPI_KEY"):
        if os.environ.get(var):
            return os.environ[var]
    rc = Path.home() / ".config" / ".pmgrc.yaml"
    if rc.exists():
        return (yaml.safe_load(rc.read_text()) or {}).get("PMG_MAPI_KEY")
    return None


def _get(params: dict) -> list[dict]:
    return _get_json(params).get("data", [])


def _get_json(params: dict) -> dict:
    key = api_key()
    if not key:
        raise MaterialsProjectError("A Materials Project API key is needed: add PMG_MAPI_KEY to ~/.config/.pmgrc.yaml.")
    try:
        r = requests.get(API, headers={"X-API-KEY": key}, params=params, timeout=60)
    except requests.RequestException as exc:
        raise MaterialsProjectError(f"Materials Project is unreachable: {type(exc).__name__}") from None
    if r.status_code != 200:
        raise MaterialsProjectError(f"Materials Project answered {r.status_code}: {r.text[:200]}")
    return r.json()


def _row(d: dict) -> dict:
    sym = d.get("symmetry") or {}
    return {"material_id": d["material_id"], "formula": d["formula_pretty"],
            "e_hull_eV_per_atom": float(d.get("energy_above_hull") or 0.0), "space_group": sym.get("symbol", ""),
            "space_group_number": sym.get("number"), "nsites": d.get("nsites"), "band_gap_eV": d.get("band_gap"),
            "theoretical": d.get("theoretical")}


def lowest_in_space_group(formula: str, space_group_number: int) -> dict:
    """The entry of `formula` in the given space group with the lowest energy above the convex hull."""
    data = _get({"formula": formula, "spacegroup_number": space_group_number, "_fields": FIELDS,
                 "_sort_fields": "energy_above_hull", "_limit": 1})
    if not data:
        raise MaterialsProjectError(f"Materials Project lists {formula} in space group {space_group_number} as empty")
    return _row(data[0])


def most_stable(formula: str) -> dict:
    data = _get({"formula": formula, "_fields": FIELDS, "_sort_fields": "energy_above_hull", "_limit": 1})
    if not data:
        raise MaterialsProjectError(f"Materials Project lists {formula} as empty")
    return _row(data[0])


def structure(material_id: str) -> dict:
    """The material's structure as Materials Project stores it (pymatgen's dict), from the cache or the API."""
    path = CACHE / f"{material_id}.json"
    if path.exists():
        return json.loads(path.read_text())
    data = _get({"material_ids": material_id, "_fields": "material_id,structure"})
    if not data:
        raise MaterialsProjectError(f"Materials Project returned an empty structure for {material_id}")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data[0]["structure"]))
    return data[0]["structure"]


def _filters(elements: list[str] | None, exclude: list[str] | None, max_elements: int | None,
             max_e_hull: float | None, max_sites: int | None) -> dict:
    """Materials Project summary filters: every listed element present, the excluded ones absent."""
    params = {}
    if elements:
        params["elements"] = ",".join(elements)
    if exclude:
        params["exclude_elements"] = ",".join(exclude)
    if max_elements:
        params["nelements_max"] = max_elements
    if max_e_hull is not None:
        params["energy_above_hull_max"] = max_e_hull
    if max_sites:
        params["nsites_max"] = max_sites
    return params


def count(elements=None, exclude=None, max_elements=None, max_e_hull=0.05, max_sites=None) -> int:
    """How many Materials Project entries match the filters."""
    meta = _get_json({**_filters(elements, exclude, max_elements, max_e_hull, max_sites), "_limit": 1,
                      "_fields": "material_id"}).get("meta", {})
    return int(meta.get("total_doc") or 0)


def search(elements=None, exclude=None, max_elements=None, max_e_hull=0.05, max_sites=None,
           limit: int | None = None, per_page: int = 500, sort: str = "material_id", formula: str | None = None) -> list[dict]:
    """Every matching entry, page by page, each row with its structure (cached as structure() caches it). `sort` is a
    Materials Project field (energy_above_hull puts the most stable first); `formula` narrows to one composition."""
    params = {**_filters(elements, exclude, max_elements, max_e_hull, max_sites), "_fields": FIELDS + ",structure",
              "_sort_fields": sort}
    if formula:
        params["formula"] = formula
    out, skip = [], 0
    per_page = per_page if limit is None else min(per_page, limit)
    while limit is None or len(out) < limit:
        data = _get({**params, "_skip": skip, "_limit": per_page})
        for d in data:
            path = CACHE / f"{d['material_id']}.json"
            if not path.exists():
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text(json.dumps(d["structure"]))
            out.append(_row(d))
        if len(data) < per_page:
            break
        skip += per_page
    return out[:limit] if limit else out

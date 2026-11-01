"""The faces of a crystal a powder can show: every stoichiometric, nonpolar, symmetric slab up to a Miller index.

A slab whose two faces match and whose composition is the bulk's has a surface energy

    γ = (E_slab − N × E_bulk per atom) / 2A

and the face with the lowest γ is the one a crystal grown near equilibrium shows most (Wulff). Slabs come from
pymatgen's SlabGenerator on the conventional standard cell, so Miller indices carry their usual meaning. A polar
slab (a dipole across it) or one whose faces differ would put a field or two different surfaces into γ, so both stay
out of the list.
"""

from __future__ import annotations

import numpy as np
from pymatgen.core import Structure
from pymatgen.core.surface import SlabGenerator, get_symmetrically_distinct_miller_indices
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer


def conventional(structure: dict) -> Structure:
    s = Structure.from_dict(structure)
    s.remove_site_property("magmom") if "magmom" in s.site_properties else None
    return SpacegroupAnalyzer(s, symprec=0.1).get_conventional_standard_structure()


def _record(slab, bulk_formula: str) -> dict:
    lat = slab.lattice.matrix
    area = float(np.linalg.norm(np.cross(lat[0], lat[1])))
    return {"hkl": [int(i) for i in slab.miller_index], "shift": round(float(slab.shift), 4),
            "n_atoms": len(slab), "area_A2": area, "formula": slab.composition.reduced_formula,
            "stoichiometric": slab.composition.reduced_formula == bulk_formula,
            "symmetric": bool(slab.is_symmetric()), "polar": bool(slab.is_polar()),
            "species": [site.specie.symbol for site in slab], "positions": slab.cart_coords.tolist(),
            "cell": lat.tolist()}


def slabs(structure: dict, max_index: int = 1, min_slab_A: float = 10.0, min_vacuum_A: float = 15.0,
          max_atoms: int = 240) -> tuple[list[dict], list[dict]]:
    """(usable, set aside): every slab of every distinct face up to `max_index`. Usable slabs are stoichiometric,
    nonpolar, symmetric and at most `max_atoms` atoms; the rest are listed with the reason in their flags."""
    bulk = conventional(structure)
    formula = bulk.composition.reduced_formula
    usable, aside = [], []
    for hkl in get_symmetrically_distinct_miller_indices(bulk, max_index):
        gen = SlabGenerator(bulk, hkl, min_slab_size=min_slab_A, min_vacuum_size=min_vacuum_A, center_slab=True,
                            primitive=True, lll_reduce=True, max_normal_search=max(hkl, key=abs) and 1)
        for slab in gen.get_slabs(symmetrize=False, ftol=0.1):
            r = _record(slab.get_orthogonal_c_slab(), formula)
            ok = r["stoichiometric"] and r["symmetric"] and not r["polar"] and r["n_atoms"] <= max_atoms
            (usable if ok else aside).append(r)
    return usable, aside


def repeats(cell: list[list[float]], width_A: float = 7.0) -> tuple[int, int]:
    """In-plane repeats that make both perpendicular widths of the cell at least `width_A`, so an adsorbate sits
    farther than that from its own periodic copies."""
    a, b = np.asarray(cell[0]), np.asarray(cell[1])
    area = np.linalg.norm(np.cross(a, b))
    wa, wb = area / np.linalg.norm(b), area / np.linalg.norm(a)
    return max(1, int(np.ceil(width_A / wa))), max(1, int(np.ceil(width_A / wb)))

"""Each ready catalyst in the measured library against the other numbers its source prints. Run:
uv run python experiments/library_consistency/check.py"""

import re
from pathlib import Path

import numpy as np
import yaml

from peroxide_engine import checks
from peroxide_engine.catalysts import library

HERE = Path(__file__).parent
data = yaml.safe_load((HERE / "data.yaml").read_text())
SUPERSCRIPT = str.maketrans("⁰¹²³⁴⁵⁶⁷⁸⁹⁻", "0123456789-")


def parse(printed: str) -> tuple[float, float]:
    """'(2.39 ± 0.09) x 10-10' -> (2.39e-10, 0.09e-10); '30' -> (30, 0); '1.1 × 107' -> (1.1e7, 0)."""
    s = printed.translate(SUPERSCRIPT).split("m2")[0].strip()
    m = re.fullmatch(r"\(?\s*([\d.]+)\s*(?:±\s*([\d.]+))?\s*\)?\s*(?:[x×]\s*10\s*(-?\d+))?", s)
    if m is None:
        raise ValueError(f"unreadable printed value {printed!r}")
    scale = 10.0 ** int(m[3]) if m[3] else 1.0
    return float(m[1]) * scale, float(m[2] or 0) * scale


def transcription(entry):
    rows = []
    for name in ("k2", "intercept", "Ea", "A", "specific_area"):
        v = getattr(entry.law, name)
        if v is None:
            continue
        value, unc = parse(v.printed)
        stated = v.uncertainty if "±" in v.printed else 0.0
        rows.append(dict(catalyst=entry.id, value_name=name, printed=v.printed, stored=v.value,
                         within=bool(np.isclose(value, v.value, rtol=1e-12) and np.isclose(unc, stated, rtol=1e-9,
                                                                                            atol=0))))
    return rows


def k1_rows(entry, spec):
    law, rows = entry.law, []
    for r in spec.get("k1_at_area", []):
        av = r["area_m2"] / (r["volume_L"] * 1e-3)
        computed = law.k_obs(av, r["T_K"])
        allowed = law.k2.uncertainty * av + law.intercept.uncertainty + r["uncertainty"]
        rows.append(dict(catalyst=entry.id, compared="k1 at a stated surface area", printed=r["value"],
                         computed=computed, allowed=allowed, within=bool(abs(computed - r["value"]) <= allowed)))
    for r in spec.get("arrhenius", []):
        ratio = law.k_obs(1.0, r["T2_K"]) / law.k_obs(1.0, r["T1_K"])
        computed = r["k1"] * ratio
        dEa = law.Ea.uncertainty * 1e3 / library.R * abs(1 / r["T1_K"] - 1 / r["T2_K"])  # in ln k
        rel = np.expm1(dEa) + r["k1_uncertainty"] / r["k1"] + r["uncertainty"] / r["value"]
        rows.append(dict(catalyst=entry.id, compared=f"k1 at {r['T2_K']} K from k1 at {r['T1_K']} K", printed=r["value"],
                         computed=computed, allowed=rel * r["value"], within=bool(abs(computed / r["value"] - 1) <= rel)))
    if "A_at_area" in spec:
        r = spec["A_at_area"]
        computed = law.A.value * np.exp(-law.Ea.value * 1e3 / (library.R * law.T_ref_K))
        rel = np.expm1(law.Ea.uncertainty * 1e3 / (library.R * law.T_ref_K)) + law.A.uncertainty / law.A.value \
            + r["k1_uncertainty"] / r["k1"]
        rows.append(dict(catalyst=entry.id, compared="A exp(-Ea/RT) at 298.15 K against k1", printed=r["k1"],
                         computed=computed, allowed=rel * r["k1"], within=bool(abs(computed / r["k1"] - 1) <= rel)))
    for r in spec.get("mass_area", []):
        computed = r["mass_g"] * law.specific_area.value
        allowed = r["mass_g"] * law.specific_area.uncertainty
        rows.append(dict(catalyst=entry.id, compared="surface area from mass and B.E.T. area", printed=r["area_m2"],
                         computed=computed, allowed=allowed, within=bool(abs(computed - r["area_m2"]) <= allowed)))
    return rows


def main():
    ready = library.runnable()
    assert {e.id for e in ready} <= set(data), "every ready catalyst needs its source's numbers in data.yaml"
    trans, consistency = [], []
    for e in ready:
        trans += transcription(e)
        consistency += k1_rows(e, data[e.id])
    passed = all(r["within"] for r in trans + consistency)
    summary = dict(
        verdict="passed" if passed else "missed",
        n_ready=len(ready), ready_names=" and ".join(e.name for e in ready),
        n_entries=len(library.entries()),
        n_awaiting=sum(e.status == "awaiting" for e in library.entries().values()),
        transcription_within=sum(r["within"] for r in trans), transcription_total=len(trans),
        consistency_within=sum(r["within"] for r in consistency), consistency_total=len(consistency),
        largest_k1_gap=max(abs(r["computed"] / r["printed"] - 1) for r in consistency if r["compared"].startswith("k1")),
        transcription=trans, consistency=consistency,
    )
    checks.write_summary(HERE, summary)
    for r in consistency:
        print(f"{r['catalyst']:<14} {r['compared']:<45} printed {r['printed']:.4g}  computed {r['computed']:.4g}  "
              f"allowed ±{r['allowed']:.2g}  {'ok' if r['within'] else 'MISS'}")
    print("verdict:", summary["verdict"])


if __name__ == "__main__":
    main()

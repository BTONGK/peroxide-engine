"""The iron model's constants and equations against the values the papers print. Run:
uv run python experiments/iron_transcription/check.py"""

from pathlib import Path

import numpy as np
import yaml

from peroxide_engine import checks, sources
from peroxide_engine.catalysts import iron

HERE = Path(__file__).parent
data = yaml.safe_load((HERE / "data.yaml").read_text())
k = iron.Constants.published()
ka = 10.0 ** -sources.value("h2o2_pKa_gallard")
h3 = iron.h_from_pH(3.0)


def round_trip(log_k: float, printed: str) -> float:
    """The authors' two roundings: log K to three decimals, then K to the printed significant figures."""
    figures = len(printed.lower().split("e")[0].replace(".", "").lstrip("0"))
    return float(f"{10.0 ** round(log_k, 3):.{figures - 1}e}")


group_a = []
for item in data["group_a"]:
    computed = iron.log_hydrolysis(item["id"], item["I_M"])
    printed = float(item["printed"])
    group_a.append(dict(name=item["name"], where=item["where"], printed=printed, computed=10.0**computed,
                        rounded_as_printed=round_trip(computed, item["printed"]),
                        log_difference=computed - np.log10(printed),
                        within=bool(round_trip(computed, item["printed"]) == printed)))

b = data["group_b"]
group_b = []


def add(name, computed, item):
    printed = float(item["printed"])
    group_b.append(dict(name=name, where=item.get("where", "De Laat & Gallard p. 2730"), printed=printed,
                        computed=computed, relative_difference=computed / printed - 1,
                        within=bool(abs(computed / printed - 1) <= 0.05)))


add("KI1 written with HO2-", k.KI1 / ka, b["KI1_prime"])
add("KI2 written with HO2-", k.KI2 / ka, b["KI2_prime"])
add("Acid constant linking the two complexes", k.KI2 * k.K1 / k.KI1, b["K_I1_I2"])
c = b["complexed_share_pH3"]
sp = iron.speciate(c["fe_M"], c["h2o2_M"], iron.h_from_pH(c["pH"]), k)
add("Share of Fe(III) in hydroperoxy complexes, pH 3, 0.5 M H2O2", (sp.I1 + sp.I2) / c["fe_M"], c)
for item in b["ki"]:
    add(f"ki at pH {item['pH']}, {item['h2o2_M']} M H2O2",
        iron.initiation_constant(iron.h_from_pH(item["pH"]), item["h2o2_M"], k), item)
k10app, k11app = iron.apparent_ho2_constants(h3, k)
add("k10app at pH 3", k10app, b["k10app_pH3"])
add("k11app at pH 3", k11app, b["k11app_pH3"])

# eq 24 as printed, and as eqs 21 and 22 combine, at pH 3 below 0.1 M H2O2
ki3 = iron.initiation_constant(h3, 0.01, k)
kd_eq21_22 = 2 * np.sqrt(ki3 * k.k7 * k11app / k10app)
kd_eq24 = np.sqrt(ki3 * k.k7 * k11app / k10app)
meas = data["kd_measured_pH3"]

passed = all(r["within"] for r in group_a) and all(r["within"] for r in group_b)
summary = dict(
    verdict="passed" if passed else "missed",
    group_a_within=sum(r["within"] for r in group_a), group_a_total=len(group_a),
    group_b_within=sum(r["within"] for r in group_b), group_b_total=len(group_b),
    group_a_max_log_difference=max(abs(r["log_difference"]) for r in group_a),
    group_b_max_relative_difference=max(abs(r["relative_difference"]) for r in group_b),
    first_run=("The first run, with the tolerance of the first amendment (0.0005 in log K), reported K1,1 at "
               "I = 0.2 M missed by 0.00007; see PASS_RULE.md."),
    kd_eq21_22=kd_eq21_22, kd_eq24=kd_eq24, kd_measured=meas["value"], kd_measured_pm=meas["plus_minus"],
    group_a=group_a, group_b=group_b,
)
checks.write_summary(HERE, summary)
for r in group_a + group_b:
    print(f"{r['name']:<58} printed {r['printed']:.4g}  computed {r['computed']:.5g}  {'ok' if r['within'] else 'MISS'}")
print(f"kd at pH 3: eqs 21-22 {kd_eq21_22:.3f}, eq 24 as printed {kd_eq24:.3f}, measured 0.47 ± 0.03")
print("verdict:", summary["verdict"])

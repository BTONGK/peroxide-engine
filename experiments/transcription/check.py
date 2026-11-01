"""Every source value as transcribed, against NIST's digitized NBS tables and each equation's printed check values.
Run: uv run python experiments/transcription/check.py (the NBS part needs cache/nbs/NBS_Tables_with_ErratumValues.csv)"""

import csv
import hashlib
import re
from pathlib import Path

from peroxide_engine import checks, oxygen, sources, water

HERE = Path(__file__).parent
CSV = checks.ROOT / "cache" / "nbs" / "NBS_Tables_with_ErratumValues.csv"
SHA = checks.ROOT / "cache" / "nbs" / "csv.sha256"
COLS = {"M": 5, "dfH": 7, "dfG": 8, "S": 10, "Cp": 11}

nbs = sources.table("nbs")
lines = list(csv.reader(CSV.open(encoding="utf-8-sig")))
clean = lambda s: re.sub(r"<[^>]+>", "", s).strip()  # noqa: E731


def csv_row(key):
    formula, rest = key.split("(", 1)
    rest = rest[:-1]
    for r in lines:
        if len(r) < 12 or clean(r[0]) != formula:
            continue
        if rest.startswith("in "):
            if r[1].strip() == rest.replace(" H2O", " H2O"):
                return r
        elif r[4].strip() == rest and not r[1].strip():
            return r
    raise KeyError(key)


compared, mismatches = 0, []
for row in nbs.rows:
    r = csv_row(row.key)
    for col, i in COLS.items():
        mine = getattr(row, col, None)
        if mine is None:
            continue
        compared += 1
        if str(mine) != r[i].strip():
            mismatches.append(dict(row=row.key, column=col, yaml=str(mine), csv=r[i].strip()))
digest = hashlib.sha256(CSV.read_bytes()).hexdigest()
published = SHA.read_text().split()[0] if SHA.exists() else None

funcs = {"iapws_vapor_pressure": water.vapor_pressure, "iapws_liquid_density": water.liquid_density,
         "o2_mole_fraction_bkp1979": oxygen.mole_fraction_1atm}
eq_rows = []
for id_, f in funcs.items():
    rec = sources.get(id_)
    for cv in rec.check_values:
        got = float(f(cv["T_K"]))
        rel = got / cv["value"] - 1
        eq_rows.append(dict(id=id_, T_K=cv["T_K"], printed=cv["value"], computed=got, rel_diff=rel,
                            tolerance=rec.check_rel_tolerance, within=bool(abs(rel) <= rec.check_rel_tolerance)))

passed = not mismatches and digest == published and all(e["within"] for e in eq_rows)
summary = dict(verdict="passed" if passed else "missed", nbs_values_compared=compared, nbs_mismatches=mismatches,
               csv_sha256=digest, csv_sha256_matches_nist=digest == published,
               equation_checks=len(eq_rows), equation_misses=sum(not e["within"] for e in eq_rows),
               max_equation_rel_diff=max(abs(e["rel_diff"]) for e in eq_rows), equations=eq_rows)
print(checks.write_summary(HERE, summary))
print(f"NBS values compared {compared}, mismatches {mismatches}, sha ok {digest == published}")
for e in eq_rows:
    print(f"{e['id']:26s} {e['T_K']:9.4f} printed {e['printed']:.6g} computed {e['computed']:.6g} rel {e['rel_diff']:+.1e} {e['within']}")
print(summary["verdict"])

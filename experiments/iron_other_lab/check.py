"""Iron(III) rates and ferrous iron measured by Barb et al. (1951), against the model as published. Run:
uv run python experiments/iron_other_lab/check.py"""

from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import yaml

from peroxide_engine import checks
from peroxide_engine.catalysts import iron

HERE = Path(__file__).parent
data = yaml.safe_load((HERE / "data.yaml").read_text())
K11A = {"used": 1.0e3, "lower": 0.0, "upper": 2.0e3}


def kd_oxygen(args):
    """Rate at nearly constant H2O2: kd where Fe(II) peaks within the first `window` of decomposition."""
    acid, fe, c0, k11a, window = args
    k = iron.Constants.published(k11a=k11a)
    h = iron.h_from_acid(acid, fe, c0, k)
    k10app, k11app = iron.apparent_ho2_constants(h, k)
    guess = 2 * np.sqrt(iron.initiation_constant(h, c0, k) * k.k7 * k11app / k10app) * fe  # steady-state kobs
    res = iron.simulate(c0, fe, h, 100 * window / guess, k, stop_at_removal=window)
    tt = np.linspace(0.0, res.sol.t[-1], 4001)
    y = res.at(tt)
    i = int(np.argmax(y[1]))
    rate = res.rate(tt[i:i + 1])[0]
    return dict(kd=rate / (y[0, i] * fe), h=h, t_peak_s=float(tt[i]))


def kd_titration(args):
    acid, fe, c0, k11a, removal = args
    k = iron.Constants.published(k11a=k11a)
    h = iron.h_from_acid(acid, fe, c0, k)
    return dict(kd=iron.kobs_initial(c0, fe, h, k, removal=removal) / fe, h=h)


def fe2(args):
    acid, fe3, c0, added, k11a = args
    k = iron.Constants.published(k11a=k11a)
    h = iron.h_from_acid(acid, fe3, c0, k)
    res = iron.simulate(c0, fe3 + added, h, 600.0, k, fe2_0_M=added)
    return dict(h=h, **{f"t{t}": float(res.at(float(t))[1, 0]) for t in (60, 180, 600)})


def score(pred, meas):
    err = np.log10(np.asarray(pred) / np.asarray(meas))
    within = float(np.mean(np.abs(err) <= 0.3))
    median = float(np.median(np.abs(err)))
    return dict(share_within_2x=within, median_abs_log10=median, max_abs_log10=float(np.max(np.abs(err))),
                passes=bool(within >= 0.8 and median <= 0.15), log10_errors=err.tolist())


def main():
    sol = data["solutions"]
    ox, ti, fx = data["kd_25C"]["oxygen_evolution"], data["kd_25C"]["titration"], data["fe2_25C"]
    out = {}
    with ProcessPoolExecutor() as pool:
        for label, k11a in K11A.items():
            kd_rows = []
            for kind, block in (("oxygen evolution", ox), ("titration", ti)):
                for s, v in block["values"].items():
                    kd_rows.append(dict(kind=kind, solution=s, acid_M=sol[s]["acid_M"], fe_M=sol[s]["fe_M"],
                                        h2o2_M=block["h2o2_M"], measured=float(v) * 1e-3))
            o_args = [(r["acid_M"], r["fe_M"], r["h2o2_M"], k11a, w) for r in kd_rows if r["kind"] == "oxygen evolution"
                      for w in (0.01, 0.005, 0.02)]
            t_args = [(r["acid_M"], r["fe_M"], r["h2o2_M"], k11a, f) for r in kd_rows if r["kind"] == "titration"
                      for f in (0.5, 0.2)]
            o_res = list(pool.map(kd_oxygen, o_args))
            t_res = list(pool.map(kd_titration, t_args))
            for r, a, b, c in zip([r for r in kd_rows if r["kind"] == "oxygen evolution"], o_res[0::3], o_res[1::3],
                                  o_res[2::3]):
                r.update(predicted=a["kd"], free_h_M=a["h"], at_half_percent=b["kd"], at_two_percent=c["kd"],
                         t_peak_s=a["t_peak_s"])
            for r, a, b in zip([r for r in kd_rows if r["kind"] == "titration"], t_res[0::2], t_res[1::2]):
                r.update(predicted=a["kd"], free_h_M=a["h"], first_20_percent=b["kd"])
            f_args = [(fx["acid_M"], row["fe_mM"] * 1e-3, row["h2o2_M"], fx["fe2_added_M"], k11a) for row in fx["rows"]]
            f_res = list(pool.map(fe2, f_args))
            fe_rows = [dict(expt=row["expt"], fe_mM=row["fe_mM"], h2o2_M=row["h2o2_M"],
                            measured=float(row["fe2_uM"]) * 1e-6, predicted=f["t180"], at_1_min=f["t60"],
                            at_10_min=f["t600"], free_h_M=f["h"]) for row, f in zip(fx["rows"], f_res)]
            kd_score = score([r["predicted"] for r in kd_rows], [r["measured"] for r in kd_rows])
            fe_score = score([r["predicted"] for r in fe_rows], [r["measured"] for r in fe_rows])
            out[label] = dict(k11a=k11a, kd=kd_score, fe2=fe_score, kd_rows=kd_rows, fe2_rows=fe_rows,
                              passes=bool(kd_score["passes"] and fe_score["passes"]))
            print(label, k11a, "kd", round(kd_score["median_abs_log10"], 3), kd_score["share_within_2x"],
                  "| Fe(II)", round(fe_score["median_abs_log10"], 3), fe_score["share_within_2x"], flush=True)
    used = out["used"]
    passed = all(v["passes"] for v in out.values())
    summary = dict(
        verdict="passed" if passed else "missed",
        kd_n=len(used["kd_rows"]), fe2_n=len(used["fe2_rows"]),
        kd_within_n=round(used["kd"]["share_within_2x"] * len(used["kd_rows"])),
        fe2_within_n=round(used["fe2"]["share_within_2x"] * len(used["fe2_rows"])),
        kd_median_factor=10 ** used["kd"]["median_abs_log10"], fe2_median_factor=10 ** used["fe2"]["median_abs_log10"],
        kd_median=used["kd"]["median_abs_log10"], kd_within=used["kd"]["share_within_2x"],
        kd_max_factor=10 ** used["kd"]["max_abs_log10"],
        fe2_median=used["fe2"]["median_abs_log10"], fe2_within=used["fe2"]["share_within_2x"],
        fe2_max_factor=10 ** used["fe2"]["max_abs_log10"],
        kd_median_lower=out["lower"]["kd"]["median_abs_log10"], kd_median_upper=out["upper"]["kd"]["median_abs_log10"],
        fe2_median_lower=out["lower"]["fe2"]["median_abs_log10"],
        fe2_median_upper=out["upper"]["fe2"]["median_abs_log10"],
        kd_rows=used["kd_rows"], fe2_rows=used["fe2_rows"],
        by_k11a=out,
    )
    checks.write_summary(HERE, summary)
    for r in used["kd_rows"]:
        print(f"{r['kind']:<17} {r['solution']}  measured {r['measured']:.4g}  model {r['predicted']:.4g}")
    for r in used["fe2_rows"]:
        print(f"Fe(II) {r['expt']:>2}  measured {r['measured']:.3g}  model {r['predicted']:.3g}")
    print("verdict:", summary["verdict"])


if __name__ == "__main__":
    main()

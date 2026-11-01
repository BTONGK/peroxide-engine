"""Iron(III) rates against De Laat & Gallard's own runs, with their one fitted constant k6 held out in turn. Run:
uv run python experiments/iron_same_lab/check.py"""

from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import yaml
from scipy.optimize import minimize_scalar

from peroxide_engine import checks
from peroxide_engine.catalysts import iron

HERE = Path(__file__).parent
data = yaml.safe_load((HERE / "data.yaml").read_text())
RUNS = ([dict(table=3, run=r["run"], pH=r["pH"], fe_M=200e-6, c0_M=10e-3, measured=float(r["kobs_1e6"]) * 1e-6)
         for r in data["table3"]]
        + [dict(table=4, run=r["run"], pH=3.0, fe_M=r["fe_uM"] * 1e-6, c0_M=r["h2o2_mM"] * 1e-3,
                measured=float(r["kobs_1e6"]) * 1e-6) for r in data["table4"]])
K11A = {"used": 1.0e3, "lower": 0.0, "upper": 2.0e3}


def kobs(args):
    run, k6, k11a = args
    k = iron.Constants.published(k6=k6, k11a=k11a)
    return iron.kobs_initial(run["c0_M"], run["fe_M"], iron.h_from_pH(run["pH"]), k)


def predict(pool, runs, k6, k11a):
    return np.array(list(pool.map(kobs, [(r, k6, k11a) for r in runs])))


def score(pred, runs):
    err = np.log10(pred / np.array([r["measured"] for r in runs]))
    within = float(np.mean(np.abs(err) <= 0.3))
    median = float(np.median(np.abs(err)))
    return dict(share_within_2x=within, median_abs_log10=median, max_abs_log10=float(np.max(np.abs(err))),
                passes=bool(within >= 0.8 and median <= 0.15), log10_errors=err.tolist())


def fit_k6(pool, runs, k11a):
    meas = np.log10([r["measured"] for r in runs])

    def loss(log_k6):
        return float(np.sum((np.log10(predict(pool, runs, 10.0**log_k6, k11a)) - meas) ** 2))

    res = minimize_scalar(loss, bounds=(-4.0, -1.0), method="bounded", options={"xatol": 2e-4})
    return float(10.0**res.x)


def main():
    t3 = [r for r in RUNS if r["table"] == 3]
    t4 = [r for r in RUNS if r["table"] == 4]
    out = {}
    with ProcessPoolExecutor() as pool:
        for label, k11a in K11A.items():
            published = predict(pool, RUNS, 2.7e-3, k11a)
            k6_t3, k6_t4, k6_all = fit_k6(pool, t3, k11a), fit_k6(pool, t4, k11a), fit_k6(pool, RUNS, k11a)
            pred_t4 = predict(pool, t4, k6_t3, k11a)
            pred_t3 = predict(pool, t3, k6_t4, k11a)
            out[label] = dict(
                k11a=k11a, published=score(published, RUNS), published_kobs=published.tolist(),
                k6_from_table3=k6_t3, table4_predicted=score(pred_t4, t4), table4_kobs=pred_t4.tolist(),
                k6_from_table4=k6_t4, table3_predicted=score(pred_t3, t3), table3_kobs=pred_t3.tolist(),
                k6_from_all=k6_all)
            r = out[label]
            r["passes"] = bool(r["published"]["passes"] and r["table4_predicted"]["passes"]
                               and r["table3_predicted"]["passes"])
            print(label, k11a, "published", r["published"]["median_abs_log10"], r["published"]["share_within_2x"],
                  "| k6 t3", k6_t3, "-> t4", r["table4_predicted"]["median_abs_log10"],
                  r["table4_predicted"]["share_within_2x"], "| k6 t4", k6_t4, "-> t3",
                  r["table3_predicted"]["median_abs_log10"], r["table3_predicted"]["share_within_2x"],
                  "| k6 all", k6_all, flush=True)
    used = out["used"]
    passed = all(v["passes"] for v in out.values())
    rows = [dict(table=r["table"], run=r["run"], pH=r["pH"], fe_uM=r["fe_M"] * 1e6, h2o2_mM=r["c0_M"] * 1e3,
                 measured_1e6=r["measured"] * 1e6, published_model_1e6=p * 1e6)
            for r, p in zip(RUNS, used["published_kobs"])]
    for row, p in zip([x for x in rows if x["table"] == 4], used["table4_kobs"]):
        row["held_out_1e6"] = p * 1e6
    for row, p in zip([x for x in rows if x["table"] == 3], used["table3_kobs"]):
        row["held_out_1e6"] = p * 1e6
    summary = dict(
        verdict="passed" if passed else "missed",
        n_runs=len(RUNS),
        published_median=used["published"]["median_abs_log10"],
        published_within=used["published"]["share_within_2x"],
        published_max_factor=10 ** used["published"]["max_abs_log10"],
        k6_from_table3=used["k6_from_table3"], k6_from_table4=used["k6_from_table4"],
        k6_from_all=used["k6_from_all"],
        t4_median=used["table4_predicted"]["median_abs_log10"], t4_within=used["table4_predicted"]["share_within_2x"],
        t4_max_factor=10 ** used["table4_predicted"]["max_abs_log10"],
        t3_median=used["table3_predicted"]["median_abs_log10"], t3_within=used["table3_predicted"]["share_within_2x"],
        t3_max_factor=10 ** used["table3_predicted"]["max_abs_log10"],
        by_k11a=out, rows=rows,
    )
    checks.write_summary(HERE, summary)
    print("verdict:", summary["verdict"])


if __name__ == "__main__":
    main()

"""Model checks: each comparison with published measurements, its rule and its verdict."""

import pandas as pd
import streamlit as st

from peroxide_engine import checks
from peroxide_engine.ui.common import verdict_badge

LABELS = {
    "rows": "Values compared",
    "leave_one_out": "Dilution series, each point left out and predicted",
    "air": "Air saturation against the USGS formula",
    "equations": "Printed check values",
    "modes": "Each gas arrangement",
    "limit": "Fast transfer against equilibrium",
    "heat": "Heat against its bookkeeping",
    "group_a": "Ionic-strength equations, as printed",
    "group_b": "Values the papers derive and print",
    "kd_rows": "Second-order constant kd (k11a = 1e3)",
    "fe2_rows": "Ferrous iron (k11a = 1e3)",
    "secondary_rows": "UMA-S 1.2.1, oxide mode (reported beside the verdict)",
}


def table(rows: list[dict]) -> pd.DataFrame:
    """Rows as a frame, floats as text to six significant figures: 7.62e-07 keeps its digits (the grid would round it
    to 0.0000007), float noise (49.99999999999999) reads as 50, and a column a row lacks stays blank."""
    return pd.DataFrame(rows).map(lambda v: ("" if v != v else f"{v:.6g}") if isinstance(v, float) else v)


class Fields(dict):
    def __missing__(self, key):
        return "{" + key + "}"


st.title("Model checks")
st.caption("Each check compares the engine with published measurements under a rule written before the comparison.")
for c in checks.all_checks():
    r, s = c["report"], c["summary"]
    with st.container(border=True):
        st.subheader(r["title"])
        st.markdown(verdict_badge(s["verdict"]))
        for label, key in (("Question", "question"), ("Compared", "compared"), ("Result", "result"),
                           ("Why it matters", "why_it_matters")):
            st.markdown(f"**{label}.** {r[key].format_map(Fields(s))}")
        tables = {k: v for k, v in s.items() if isinstance(v, list) and v and isinstance(v[0], dict)}
        with st.expander("Numbers"):
            for key, rows in tables.items():
                st.markdown(f"**{LABELS.get(key, key)}**")
                st.dataframe(table(rows), hide_index=True, width="stretch")
        with st.expander("Pass rule"):
            st.markdown((checks.EXPERIMENTS / c["folder"] / "PASS_RULE.md").read_text())

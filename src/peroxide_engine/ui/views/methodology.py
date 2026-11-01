"""Methodology: how each number is computed, every source value with its page, the limits and the audit."""

import pandas as pd
import streamlit as st

from peroxide_engine import sources
from peroxide_engine.ui.common import doc

m = doc("methodology")
st.title("Methodology")
st.markdown(m["intro"])
tabs = st.tabs(["How it works", "Every source value", "What comes later", "Audit"])
with tabs[0]:
    for sec in m["sections"]:
        st.subheader(sec["title"])
        st.markdown(sec["body"])
        for eq in sec.get("equations", []):
            st.code(eq, language=None)
        if "uncertainty" in sec:
            st.markdown(f"**Uncertainty.** {sec['uncertainty']}")
with tabs[1]:
    reg = sources.registry()
    consts = [r for r in reg.values() if r.kind == "constant"]
    st.subheader("Constants")
    st.dataframe(pd.DataFrame([{"Name": r.name, "Value": f"{r.value:.15g}", "Units": r.units,
                                "Uncertainty": "exact" if r.exact else f"{r.uncertainty:.15g}",
                                "Source": r.source.cite,
                                "Where": r.source.locator} for r in consts]), hide_index=True, width="stretch")
    for t in (r for r in reg.values() if r.kind == "table"):
        st.subheader(t.name)
        st.caption(t.uncertainty_note)
        cols = [c for c in ("M", "dfH", "dfG", "S", "Cp", "n_water")]
        st.dataframe(pd.DataFrame([{"Substance": row.key, **{c: getattr(row, c, "") for c in cols},
                                    "Where": row.locator} for row in t.rows]), hide_index=True, width="stretch")
        st.caption(f"{t.source.reference}. Columns: " + "; ".join(f"{k}: {v}" for k, v in t.columns.items()))
    st.subheader("Equations")
    for c in (r for r in reg.values() if r.kind == "correlation"):
        with st.container(border=True):
            st.markdown(f"**{c.name}** ({c.source.cite}, {c.source.locator})")
            st.code(c.equation, language=None)
            st.markdown("Coefficients: " + ", ".join(f"{k} = {v:.12g}" for k, v in c.coefficients.items()))
            st.caption(f"Accuracy: {c.accuracy_note} Read from: {c.transcribed_from}.")
    st.subheader("References")
    seen = {}
    for r in reg.values():
        seen.setdefault(r.source.reference, r.source)
    for ref, src in seen.items():
        link = f"https://doi.org/{src.doi}" if src.doi else src.url
        st.markdown(f"- {ref}. [{src.doi or 'link'}]({link})")
with tabs[2]:
    st.markdown(m["limits_intro"])
    for item in m["limits"]:
        st.markdown(f"**{item['title']}.** {item['body']}")
with tabs[3]:
    for a in m["audit"]:
        st.markdown(f"- {a}")

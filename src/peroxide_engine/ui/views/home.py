"""Home page."""

import streamlit as st

from peroxide_engine import checks, sources
from peroxide_engine.ui import pages
from peroxide_engine.ui.common import doc

st.title("Peroxide engine")
st.markdown(
    "Computes what an experimentalist measures when hydrogen peroxide decomposes: the rate, the oxygen in the liquid "
    "and the gas, the pressure, the heat and the titration volumes. Every constant traces to a page of a published "
    "source, and every model reaches this app only after it reproduces measurements under a rule written down first."
)
results = checks.all_checks()
passed = sum(c["summary"]["verdict"] == "passed" for c in results)
reg = sources.registry()
n_values = sum(len(r.rows) if r.kind == "table" else 1 for r in reg.values())
m1, m2, m3 = st.columns(3)
m1.metric("Model checks passed", f"{passed} of {len(results)}")
m2.metric("Source values, each with its page", n_values)
m3.metric("Phase", "1 of 6")

st.subheader("What you can compute today")
st.markdown(
    "- **Any catalyst whose first-order rate constant you know:** the full time course of H2O2, oxygen formed, "
    "dissolved oxygen, gas burette volume or sealed-vessel pressure, and titration volumes.\n"
    "- **Heat and heat flow** at 25 °C, as an isothermal calorimeter records them, including the heat taken as oxygen "
    "leaves the liquid.\n"
    "- **The lab's fits:** rate constant, half-life and initial rate read from the computed curves."
)
st.page_link(pages.page("run"))

st.subheader("What comes next")
for ph in doc("status")["phases"]:
    if ph["state"] != "done":
        st.markdown(f"- **Phase {ph['n']}, {ph['name']}:** {ph['summary']}")
st.page_link(pages.page("status"))

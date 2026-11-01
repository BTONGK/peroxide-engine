"""Catalyst library: every catalyst the engine runs or will run, by class, with its source and measured ranges."""

import pandas as pd
import streamlit as st

from peroxide_engine.catalysts import library

st.title("Catalyst library")
st.caption("Each catalyst runs on the Run page inside the conditions its laboratory measured. Catalysts waiting for a "
           "paper are listed with the paper they need.")

entries = library.entries().values()
ready = [e for e in entries if e.status == "ready"]
c1, c2, c3 = st.columns(3)
c1.metric("Ready to run", len(ready))
c2.metric("Waiting for a paper", len(entries) - len(ready))
c3.metric("Catalyst classes", len({e.cls for e in entries}))


def measured(e: library.Entry) -> str:
    law = e.law
    if isinstance(law, library.SurfaceFirstOrder):
        T = "–".join(f"{t - 273.15:.0f}" for t in law.T_range_K)
        c = "–".join(dict.fromkeys(f"{x * 1e3:g}" for x in law.c0_M))
        a = "–".join(f"{x / 1e3:,.0f}" for x in law.area_per_volume_per_m)
        return f"{T} °C; H2O2 {c} mM; {a} m² of oxide per liter"
    if isinstance(law, library.Mechanism):
        return law.summary
    return ""


for cls, title in library.CLASSES.items():
    group = [e for e in entries if e.cls == cls]
    st.subheader(title)
    rows = []
    for e in group:
        law = e.law
        rate = (f"k2 = {law.k2.printed} m/s, Ea = {law.Ea.printed} kJ/mol"
                if isinstance(law, library.SurfaceFirstOrder) else "")
        rows.append({"Catalyst": e.name, "Status": "Ready" if e.status == "ready" else "Waiting for a paper",
                     "Source": e.source.cite, "Link": f"https://doi.org/{e.source.doi}" if e.source.doi else e.source.url,
                     "Measured rate": rate, "Measured over": measured(e), "Needs": e.needs or "",
                     "Material": e.material})
    st.dataframe(pd.DataFrame(rows), hide_index=True, width="stretch",
                 column_config={"Link": st.column_config.LinkColumn("Link", display_text="open")})
    for e in group:
        if e.status == "ready" and isinstance(e.law, library.SurfaceFirstOrder):
            with st.expander(f"{e.name}: where each number comes from"):
                law = e.law
                for name in ("k2", "intercept", "Ea", "A", "specific_area"):
                    v = getattr(law, name)
                    if v is not None:
                        st.markdown(f"**{name}** = {v.printed} {v.units} — {v.locator}" + (f". {v.note}" if v.note else ""))
                st.markdown(f"**Ranges.** {law.ranges_from}")
                if e.notes:
                    st.markdown(e.notes)

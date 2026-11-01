"""Project status, from docs/status.yaml."""

import pandas as pd
import streamlit as st

from peroxide_engine.ui.common import doc

s = doc("status")
st.title("Project status")
st.caption(f"Updated {s['updated']}")
st.markdown(s["summary"])
BADGE = {"done": ":green-badge[Done]", "next": ":orange-badge[Next]", "planned": ":gray-badge[Planned]"}
st.subheader("Phases")
for ph in s["phases"]:
    st.markdown(f"{BADGE[ph['state']]} **Phase {ph['n']}, {ph['name']}.** {ph['summary']}")
st.subheader("Recent work")
for r in s["recent"]:
    st.markdown(f"- **{r['date']}** {r['text']}")
st.subheader("Next")
for n in s["next"]:
    st.markdown(f"- {n}")
st.subheader("Papers to fetch")
st.markdown(s["papers"]["intro"])
st.dataframe(pd.DataFrame([{"Phase": p["phase"], "Paper": p["cite"], "DOI": f"https://doi.org/{p['doi']}",
                            "What the engine takes from it": p["need"]} for p in s["papers"]["list"]]),
             hide_index=True, width="stretch", column_config={"DOI": st.column_config.LinkColumn("DOI")})
st.subheader("Findings")
for f in s["findings"]:
    st.markdown(f"- {f}")
st.subheader("Decisions")
for d in s["decisions"]:
    st.markdown(f"- {d}")

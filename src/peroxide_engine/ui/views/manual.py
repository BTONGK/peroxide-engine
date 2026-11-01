"""Manual."""

import streamlit as st

from peroxide_engine.ui.common import doc

m = doc("manual")
st.title("Manual")
st.markdown(m["intro"])
for step in m["steps"]:
    st.subheader(step["title"])
    st.markdown(step["body"])
st.subheader("Terms")
for g in m["glossary"]:
    st.markdown(f"**{g['term']}:** {g['meaning']}")

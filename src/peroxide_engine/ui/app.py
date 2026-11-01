"""Peroxide engine app. Launch with `uv run pe ui`."""

import streamlit as st

from peroxide_engine.ui import pages

st.set_page_config(page_title="Peroxide engine", page_icon=":material/science:", layout="wide")
st.navigation(pages.sections()).run()

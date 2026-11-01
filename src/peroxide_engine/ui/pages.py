"""The app's pages, shared by the navigation and the links between pages."""

from pathlib import Path

import streamlit as st

VIEWS = Path(__file__).parent / "views"


def page(name: str) -> st.Page:
    title, icon = TITLES[name]
    return st.Page(VIEWS / f"{name}.py", title=title, icon=icon, default=name == "home")


TITLES = {
    "home": ("Home", ":material/home:"),
    "manual": ("Manual", ":material/menu_book:"),
    "run": ("Run an experiment", ":material/experiment:"),
    "checks": ("Model checks", ":material/fact_check:"),
    "methodology": ("Methodology", ":material/library_books:"),
    "status": ("Project status", ":material/dashboard:"),
}


def sections() -> dict[str, list[st.Page]]:
    return {
        "Start": [page("home"), page("manual")],
        "Experiment": [page("run")],
        "Accuracy": [page("checks"), page("methodology")],
        "Project": [page("status")],
    }

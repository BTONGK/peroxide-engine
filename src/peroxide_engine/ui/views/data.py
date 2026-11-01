"""Data collection: every computed crystal, face and binding energy, and the measured values they are tested against."""

import io
import json

import pandas as pd
import streamlit as st

from peroxide_engine.compute import collection


@st.cache_data
def load(_fingerprint: tuple) -> dict[str, pd.DataFrame]:
    return collection.tables()


def narrow(df: pd.DataFrame, text: str, models: list[str], origins: list[str]) -> pd.DataFrame:
    if df.empty:
        return df
    keep = df["Model"].isin(models) & df["Origin"].isin(origins)
    if text:
        keep &= df["Formula"].str.contains(text, case=False, regex=False) | df["Materials Project"].str.contains(
            text, case=False, regex=False)
    return df[keep]


def csv(df: pd.DataFrame) -> bytes:
    return df.to_csv(index=False).encode()


st.title("Data collection")
st.caption("Every crystal the engine has computed, from runs on the Compute a catalyst page, screens prepared with "
           "pe screen and the model checks, kept in one place with the model and run that made each number. The "
           "Measured tab holds the laboratory values the computed data are tested against.")

data = load(collection.fingerprint())
crystals, faces, sites = data["crystals"], data["faces"], data["sites"]
c1, c2, c3, c4 = st.columns(4)
c1.metric("Crystals computed", crystals["Materials Project"].nunique() if len(crystals) else 0)
c2.metric("Faces relaxed", len(faces))
c3.metric("Binding energies", len(sites))
c4.metric("Runs", crystals["Run"].nunique() if len(crystals) else 0)

if len(crystals):
    f1, f2, f3 = st.columns([2, 2, 2])
    text = f1.text_input("Formula or Materials Project ID contains", "")
    models = f2.multiselect("Model", sorted(crystals["Model"].unique()), default=sorted(crystals["Model"].unique()))
    origins = f3.multiselect("From", sorted(crystals["Origin"].unique()), default=sorted(crystals["Origin"].unique()))
    crystals, faces, sites = (narrow(d, text, models, origins) for d in (crystals, faces, sites))

tab_c, tab_f, tab_s, tab_m = st.tabs(["Crystals", "Faces", "Binding sites", "Measured"])
fmt = {"Surface energy (J/m²)": st.column_config.NumberColumn(format="%.3f")}
with tab_c:
    if crystals.empty:
        st.info("Crystals appear here as runs on the Compute a catalyst page finish their first face.")
    else:
        st.caption("One row per crystal per run: the face of lowest surface energy and the strongest site of each "
                   "molecule on it, where the molecule stays whole. Binding energies are against gas molecules "
                   "(OH: H2O − ½ H2; O: H2O − H2; OOH: 2 H2O − 3/2 H2; H: ½ H2); negative binds.")
        st.dataframe(crystals, hide_index=True, width="stretch", column_config=fmt)
        st.download_button("Crystals (CSV)", csv(crystals), "crystals.csv", icon=":material/download:")
with tab_f:
    if not faces.empty:
        st.caption("Every stoichiometric, nonpolar, symmetric face relaxed, with its surface energy.")
        st.dataframe(faces, hide_index=True, width="stretch", column_config=fmt)
        st.download_button("Faces (CSV)", csv(faces), "faces.csv", icon=":material/download:")
with tab_s:
    if not sites.empty:
        st.caption("Every molecule on every site, as relaxed. Kept marks the sites where the molecule stayed whole "
                   "and on the surface; Flag names what happened at the others.")
        st.dataframe(sites, hide_index=True, width="stretch")
        st.download_button("Binding sites (CSV)", csv(sites), "binding_sites.csv", icon=":material/download:")
with tab_m:
    measured = collection.measured()
    st.caption("Measured values, each with its source and where the engine uses it.")
    st.dataframe(measured, hide_index=True, width="stretch")
    st.download_button("Measured (CSV)", csv(measured), "measured.csv", icon=":material/download:")

everything = {k: v.to_dict(orient="records") for k, v in
              (("crystals", crystals), ("faces", faces), ("sites", sites), ("measured", collection.measured()))}
st.download_button("Everything shown (JSON)", io.BytesIO(json.dumps(everything, indent=1, default=str).encode()),
                   "peroxide_engine_data.json", icon=":material/database:")

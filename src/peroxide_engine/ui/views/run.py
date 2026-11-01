"""Run an experiment: a catalyst's rate constant and the apparatus in, what the instruments read out."""

import altair as alt
import numpy as np
import pandas as pd
import streamlit as st

from peroxide_engine import assays, vessel
from peroxide_engine.ui.common import sig

T0 = 273.15
ARRANGEMENTS = {"Gas burette": "burette", "Sealed vessel with a pressure gauge": "closed", "Open vessel": "open"}

st.title("Run an experiment")
st.caption("Give your catalyst's rate constant and the apparatus; the engine computes what each instrument reads.")

c1, c2, c3 = st.columns(3)
with c1:
    st.subheader("Catalyst and solution")
    k = st.number_input("Rate constant k at the run temperature (1/s)", min_value=1e-9, max_value=10.0, value=1e-3,
                        format="%.3e", help="First order in H2O2. From a half-life: k = 0.693 / half-life.")
    c0_mM = st.number_input("Starting H2O2 (mmol/L)", min_value=0.001, max_value=1000.0, value=100.0, format="%.3f")
    t_C = st.number_input("Temperature (°C)", min_value=5.0, max_value=55.0, value=25.0, step=0.5, format="%.2f")
with c2:
    st.subheader("Apparatus")
    arrangement = st.radio("Gas arrangement", list(ARRANGEMENTS))
    liquid = st.number_input("Liquid volume (mL)", min_value=0.1, max_value=1e5, value=50.0)
    gas = st.number_input("Gas space above the liquid, with tubing (mL)", min_value=0.1, max_value=1e5, value=50.0)
    p_kPa = st.number_input("Room pressure, or starting pressure if sealed (kPa)", min_value=50.0, max_value=200.0,
                            value=101.325, format="%.3f")
with c3:
    st.subheader("Oxygen transfer and run length")
    transfer = st.radio("Dissolved oxygen", ["In balance with the gas (well stirred)", "Crosses the surface at kLa"])
    kla = None
    if transfer.startswith("Crosses"):
        kla = st.number_input("kLa (1/s)", min_value=1e-6, max_value=100.0, value=0.01, format="%.4f",
                              help="Your vessel's oxygen transfer coefficient, measured by the dynamic method.")
    halves = st.slider("Run length (half-lives)", 1.0, 10.0, 5.0, 0.5)

setup = vessel.Setup(c0_M=c0_mM / 1000, k_per_s=k, liquid_mL=liquid, gas_mL=gas, T_K=t_C + T0,
                     pressure_Pa=p_kPa * 1000, headspace=ARRANGEMENTS[arrangement], kLa_per_s=kla,
                     t_end_s=halves * np.log(2) / k)
problems = setup.validate()
if problems:
    for p in problems:
        st.error(p)
    st.stop()


@st.cache_data(show_spinner="Computing the run…")
def compute(s: vessel.Setup):
    run = vessel.simulate(s)
    return run.series, run.readouts, run.notes


series, readouts, notes = compute(setup)
for n in notes:
    st.warning(n)

t_end = series["time_s"].iloc[-1]
unit, div = ("h", 3600) if t_end > 7200 else ("min", 60) if t_end > 120 else ("s", 1)
df = series.copy()
df["time"] = df["time_s"] / div

st.subheader("Readouts")
table = pd.DataFrame([{"Readout": r["name"], "Value": sig(float(r["value"])), "Units": r["units"],
                       "How it is computed": r["how"]} for r in readouts])
st.dataframe(table, hide_index=True, width="stretch")
if "heat_flow_mW" not in series:
    st.caption("Heat appears for runs at 25.00 °C in the burette and sealed arrangements, where the tables measure it "
               "and the energy accounting is complete; the Methodology page explains.")


def chart(col: str, title: str, y: str):
    return (alt.Chart(df).mark_line(color="#0f766e")
            .encode(x=alt.X("time", title=f"Time ({unit})"), y=alt.Y(col, title=y))
            .properties(title=title, height=260))


gas_col = next(c for c in ("gas_collected_mL", "pressure_kPa", "O2_released_mmol") if c in series)
gas_title = {"gas_collected_mL": ("Gas collected", "mL"), "pressure_kPa": ("Pressure", "kPa"),
             "O2_released_mmol": ("O2 released to the room", "mmol")}[gas_col]
g1, g2 = st.columns(2)
g1.altair_chart(chart("H2O2_mM", "H2O2 left", "mmol/L"), width="stretch")
g2.altair_chart(chart(gas_col, gas_title[0], gas_title[1]), width="stretch")
g3, g4 = st.columns(2)
g3.altair_chart(chart("dissolved_O2_mg_per_L", "Dissolved oxygen", "mg/L"), width="stretch")
if "heat_flow_mW" in series:
    g4.altair_chart(chart("heat_flow_mW", "Heat flow (isothermal calorimeter)", "mW"), width="stretch")

with st.expander("Titrations of an aliquot", icon=":material/science:"):
    a1, a2, a3 = st.columns(3)
    aliquot = a1.number_input("Aliquot (mL)", min_value=0.01, value=1.00, format="%.2f")
    c_mn = a2.number_input("KMnO4 (mol/L)", min_value=1e-5, value=0.0200, format="%.4f")
    c_th = a3.number_input("Na2S2O3 (mol/L)", min_value=1e-5, value=0.1000, format="%.4f")
    half = np.log(2) / k
    times = [0.0, half, 2 * half, t_end]
    rows = []
    for tt in times:
        c = float(np.interp(tt, series["time_s"], series["H2O2_mM"])) / 1000
        rows.append({f"Time ({unit})": sig(tt / div), "H2O2 (mmol/L)": sig(c * 1000),
                     "KMnO4 (mL)": sig(assays.permanganate_titre_mL(c, aliquot, c_mn)),
                     "Na2S2O3 (mL)": sig(assays.thiosulfate_titre_mL(c, aliquot, c_th))})
    st.dataframe(pd.DataFrame(rows), hide_index=True, width="stretch")
    st.caption("Permanganate: 2 MnO4⁻ + 5 H2O2 + 6 H⁺ → 2 Mn²⁺ + 5 O2 + 8 H2O. Iodometric: H2O2 + 2 I⁻ + 2 H⁺ → I2 + "
               "2 H2O, then I2 + 2 S2O3²⁻ → 2 I⁻ + S4O6²⁻. Each titre assumes the aliquot is quenched the moment it "
               "is drawn.")

st.download_button("Download the time series (CSV)", series.to_csv(index=False).encode(), "peroxide_run.csv",
                   "text/csv", icon=":material/download:")

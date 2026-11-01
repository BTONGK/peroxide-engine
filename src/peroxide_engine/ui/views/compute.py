"""Compute a catalyst: pick crystals from the Materials Project and compute them with Meta's UMA model, on demand."""

import time

import pandas as pd
import streamlit as st
from ase.data import chemical_symbols

from peroxide_engine import checks
from peroxide_engine.compute import materials, screens
from peroxide_engine.ui import pages
from peroxide_engine.ui.common import verdict_badge

ELEMENTS = chemical_symbols[1:95]
REFERENCES = {"OH": "H2O − ½ H2", "O": "H2O − H2", "OOH": "2 H2O − 3/2 H2", "H": "½ H2"}
STATES = {"queued": ":gray-badge[In the queue]", "running": ":blue-badge[Computing]",
          "finished": ":green-badge[Finished]", "stopped": ":orange-badge[Stopped]", "failed": ":red-badge[Failed]"}


def hkl(v) -> str:
    return "(" + " ".join(str(i).replace("-", "−") for i in v) + ")" if v else ""


def span(seconds: tuple[float, float], about: bool = True) -> str:
    lo, hi = (s / 3600 for s in seconds)
    text = (f"{max(1, round(lo * 60))}–{max(1, round(hi * 60))} min" if hi < 1 else f"{lo:.1f}–{hi:.1f} h")
    return f"about {text}" if about else text


def elapsed(seconds: float) -> str:
    return f"{seconds / 60:.0f} min" if seconds < 3600 else f"{seconds / 3600:.1f} h"


st.title("Compute a catalyst")
st.caption("Pick crystals from the Materials Project and the engine computes them with Meta's UMA model: the faces "
           "each crystal shows, the surface energy of each face, and how strongly OH, O, OOH or H bind on every site "
           "of the most stable face. Runs compute one after another on this computer and keep going when the page "
           "closes; a cluster runs the same job folders side by side.")

check = checks.load(checks.EXPERIMENTS / "oxide_descriptor_check")["summary"]
rho = f"{check['rho_Ea']:.2f}".replace("-", "−")
with st.container(border=True):
    st.markdown(
        f"**What the numbers are.** Energies from UMA, which places OH on oxides within about 0.2 eV of DFT (OC22). "
        f"Faces and surface energies match the literature for seven of the eight oxides tested. How these energies "
        f"relate to activity is under test: OH binding on the dry face ranked eight measured oxides at "
        f"ρ = {rho} against their barriers, where the rule asks for +0.6, "
        f"{verdict_badge(check['verdict'])}. So this page reports energies, and rates come from the Run page's "
        f"checked models.")
    st.page_link(pages.page("checks"), label="The check, with every face and site")

new, past = st.tabs(["New run", "Runs"])

with new:
    how = st.radio("Find crystals", ["By formula", "By elements"], horizontal=True)
    try:
        if how == "By formula":
            c1, c2 = st.columns([2, 1])
            formula = c1.text_input("Formula", "MnO2", help="As the Materials Project writes it: MnO2, Fe2O3, LaCoO3")
            hull = c2.number_input("Highest energy above the convex hull (eV/atom)", 0.0, 1.0, 0.1, 0.01)
            if st.button("Search", key="search_formula", icon=":material/search:"):
                st.session_state["found"] = materials.search(formula=formula.strip(), max_e_hull=hull, limit=50,
                                                             sort="energy_above_hull")
                st.session_state["found_total"] = len(st.session_state["found"])
        else:
            c1, c2 = st.columns(2)
            elements = c1.multiselect("Elements every crystal contains", ELEMENTS, default=["O"])
            exclude = c2.multiselect("Elements to leave out", ELEMENTS)
            c3, c4, c5 = st.columns(3)
            most = c3.number_input("Most elements per crystal", 1, 6, 2)
            sites = c4.number_input("Most atoms in the crystal cell", 2, 200, 20)
            hull = c5.number_input("Highest energy above the convex hull (eV/atom)", 0.0, 1.0, 0.02, 0.01)
            if st.button("Search", key="search_elements", icon=":material/search:"):
                filters = dict(elements=elements, exclude=exclude, max_elements=most, max_e_hull=hull, max_sites=sites)
                st.session_state["found_total"] = materials.count(**filters)
                st.session_state["found"] = materials.search(**filters, limit=500, sort="energy_above_hull")
    except materials.MaterialsProjectError as exc:
        st.error(str(exc))

    found = st.session_state.get("found")
    if found is not None:
        total = st.session_state.get("found_total", len(found))
        st.caption(f"{len(found):,} crystals, most stable first"
                   + (f", the first {len(found):,} of {total:,}" if total > len(found) else "")
                   + ". Tick the ones to compute.")
        table = pd.DataFrame([{
            "Formula": r["formula"], "Materials Project": r["material_id"], "Space group": r["space_group"],
            "Above the hull (eV/atom)": round(r["e_hull_eV_per_atom"], 3), "Atoms in cell": r["nsites"],
            "Band gap (eV)": None if r["band_gap_eV"] is None else round(r["band_gap_eV"], 2),
            "Made in a lab": "predicted" if r["theoretical"] else "yes"} for r in found])
        event = st.dataframe(table, on_select="rerun", selection_mode="multi-row", hide_index=True, width="stretch",
                             key=f"pick-{len(found)}-{found[0]['material_id'] if found else ''}")
        picked = [found[i] for i in event.selection.rows]
        c1, c2 = st.columns([1, 2])
        adsorbates = c1.multiselect("Molecules to bind", screens.ADSORBATES, default=["OH"],
                                    help="Each is set against gas molecules: " +
                                         "; ".join(f"{a} against {r}" for a, r in REFERENCES.items()))
        default = ", ".join(r["formula"] for r in picked[:4]) + (" and more" if len(picked) > 4 else "")
        name = c2.text_input("Name of the run", default or "My run")
        if picked and adsorbates and st.button(f"Compute {len(picked)} crystal{'s' * (len(picked) > 1)}",
                                               type="primary", icon=":material/play_arrow:"):
            with st.spinner("Cutting every stable face of each crystal"):
                ready, aside = screens.prepare(picked)
            for a in aside:
                st.warning(f"{a['formula']} ({a['id']}) is set aside: {a['reason']}.")
            if ready:
                folder = screens.create(name, ready, aside, adsorbates)
                screens.start_runner()
                meta = screens.runs()[0]["meta"] if screens.runs() else {}
                st.success(f"Queued {len(ready)} crystal{'s' * (len(ready) > 1)} with "
                           f"{sum(m['faces'] for m in meta.get('materials', []))} faces, "
                           f"{span(meta.get('estimate_s', (0, 0)))} on this computer. Follow it on the Runs tab.")


def show(folder):
    run = next(r for r in screens.runs() if r["folder"] == folder)
    meta, state = run["meta"], screens.state(run)
    done, planned = screens.progress(run)
    c1, c2, c3, c4 = st.columns(4)
    c1.metric("Crystals", len(meta["materials"]))
    c2.metric("Relaxations finished", f"{done} of ~{planned}")
    c3.metric("Estimate here", span(meta["estimate_s"], about=False), help="On this computer's CPU; a GPU runs faster")
    started, updated = run["status"].get("started"), run["status"].get("updated")
    c4.metric("Computing for", elapsed(updated - started) if started and updated else "–")
    message = run["status"].get("message", "") if state in ("running", "failed") else ""
    st.markdown(f"{STATES.get(state, state)} {message}")
    st.progress(1.0 if state == "finished" else min(done / planned, 1.0))
    b1, b2, _ = st.columns([1, 1, 4])
    if state in ("queued", "running") and b1.button("Stop", icon=":material/stop:", key=f"stop-{folder.name}"):
        screens.stop(folder)
        st.rerun()
    if state in ("stopped", "failed") and b2.button("Resume", icon=":material/play_arrow:", key=f"go-{folder.name}"):
        screens.resume(folder)
        st.rerun()

    summary = screens.summary(folder)
    ads = meta["adsorbates"]
    rows = []
    for m in summary:
        row = {"Crystal": m["label"], "Materials Project": m["id"],
               "State": {"finished": "finished", "computing": "computing", "waiting": "waiting", "error": "error"}[m["state"]],
               "Face": hkl(m["face"]),
               "Surface energy (J/m²)": min((f["gamma_J_m2"] for f in m["faces"]), default=None)}
        for a in ads:
            b = m["best"].get(a)
            row[f"{a} binding (eV)"] = None if b is None else round(b["E_bind_eV"], 3)
            row[f"{a} site"] = "" if b is None else b["site"]
        if "OH" in ads:
            row["H2O2 → 2 OH* (eV)"] = None if m["D_eV"] is None else round(m["D_eV"], 3)
        rows.append(row)
    st.dataframe(pd.DataFrame(rows), hide_index=True, width="stretch",
                 column_config={"Surface energy (J/m²)": st.column_config.NumberColumn(format="%.2f")})
    st.caption("Binding energies: negative binds; each molecule against its gas reference (" +
               "; ".join(f"{a}: {REFERENCES[a]}" for a in ads) + "), at the strongest site where it stays whole. "
               "H2O2 → 2 OH*: the energy of splitting one H2O2 into two OH on the face, from the strongest OH site.")

    for m in summary:
        if not m["faces"]:
            continue
        with st.expander(f"{m['label']} · {m['id']}"):
            if m["error"]:
                st.error(m["error"])
            st.markdown("**Faces**, each relaxed by UMA; the lowest surface energy is the face used")
            st.dataframe(pd.DataFrame([{"Face": hkl(f["hkl"]), "Termination shift": f["shift"], "Atoms": f["atoms"],
                                        "Surface energy (J/m²)": round(f["gamma_J_m2"], 3),
                                        "Used": "yes" if f["index"] == m["face_index"] else "",
                                        "Converged": "yes" if f["converged"] else "at the step limit"}
                                       for f in m["faces"]]), hide_index=True, width="stretch")
            if m["sites"]:
                st.markdown(f"**Sites** on {hkl(m['face'])}, repeated {m['repeat'][0]}×{m['repeat'][1]} "
                            f"({m['face_atoms']} atoms)")
                st.dataframe(pd.DataFrame([{"Molecule": s["adsorbate"], "Site": s["site"],
                                            "Binding (eV)": round(s["E_bind_eV"], 3),
                                            "Kept": "yes" if s["kept"] else s["anomaly"],
                                            "Converged": "yes" if s["converged"] else "at the step limit"}
                                           for s in m["sites"]]), hide_index=True, width="stretch")
            files = [("Face", f"structures/{m['id'].replace(':', '_').replace('/', '_')}__face.cif")]
            files += [(f"{a} at {b['site']}", next(s["structure"] for s in m["sites"]
                                                   if s["adsorbate"] == a and s["site"] == b["site"]))
                      for a, b in m["best"].items()]
            cols = st.columns(len(files))
            for col, (label, rel) in zip(cols, files):
                path = folder / rel
                if path.exists():
                    col.download_button(f"{label} (CIF)", path.read_bytes(), file_name=f"{m['formula']}_{path.name}",
                                        key=f"dl-{folder.name}-{rel}")
    if meta["set_aside"]:
        st.caption("Set aside: " + "; ".join(f"{a['formula']} ({a['id']}): {a['reason']}" for a in meta["set_aside"]))


with past:
    every = screens.runs()
    if not every:
        st.info("Runs you start on the New run tab appear here, with their progress and results.")
    else:
        names = {f"{r['meta']['name']} · {time.strftime('%d %b %H:%M', time.localtime(r['meta']['created']))}": r
                 for r in every}
        pick = st.selectbox("Run", list(names))
        chosen = names[pick]
        live = screens.state(chosen) in ("queued", "running")
        st.fragment(show, run_every=10 if live else None)(chosen["folder"])

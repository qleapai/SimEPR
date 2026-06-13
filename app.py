from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd
import streamlit as st

from epr_simfit.about import ABOUT_TEXT, DETAILED_INFO, SHORT_CREDIT
from epr_simfit.constants import DEFAULT_MW_FREQUENCY_GHZ
from epr_simfit.demo_data import DEMO_CONDITIONS, generate_demo_text
from epr_simfit.export import build_export_zip, fit_components_dataframe
from epr_simfit.fitter import fit_spectrum
from epr_simfit.io import parse_epr_text
from epr_simfit.metadata_parser import metadata_table
from epr_simfit.model_comparison import compare_models
from epr_simfit.model_library import MODEL_DESCRIPTIONS, MODEL_PRESETS, component_table, default_components
from epr_simfit.model_suggester import ExperimentContext, suggest_models
from epr_simfit.plotting import comparison_bar_figure, component_figure, residual_figure, spectrum_figure
from epr_simfit.interpretation import (
    SCIENCE_FIT,
    SCIENCE_IMPORT,
    SCIENCE_MODELBUILDER,
    SCIENCE_PREPROCESS,
    assess_fit_quality,
    intermediates_dataframe,
    publication_methods_paragraph,
    publication_parameters_table,
    suggest_fit_improvements,
    suggest_intermediates,
)
from epr_simfit.preprocessing import preprocess_spectrum
from epr_simfit.report import generate_report_html, generate_report_text
from epr_simfit.simulator import simulate_model
from epr_simfit.spin_models import Nucleus, SpinComponent


APP_DIR = Path(__file__).resolve().parent
WHITE_PAPER = APP_DIR / "docs" / "WHITE_PAPER.md"
CITATION = APP_DIR / "CITATION.cff"
LOGO     = APP_DIR / "assets" / "logo.svg"
LOGO_BANNER = APP_DIR / "assets" / "logo_banner.svg"

st.set_page_config(page_title="SimEPR", page_icon="⚛️", layout="wide")

st.markdown(
    """
    <style>
    .block-container {padding-top: 1.35rem; padding-bottom: 2rem;}
    [data-testid="stMetric"] {background:#f8fafc;border:1px solid #e5e7eb;border-radius:8px;padding:0.65rem 0.8rem;}
    .credit-note {font-size:0.9rem;color:#44515c;margin-top:-0.35rem;margin-bottom:1rem;}
    .science-note {background:#eef6ff;border:1px solid #bfdbfe;border-radius:8px;padding:0.8rem 1rem;color:#1e3a5f;}
    .warning-note {background:#fff7ed;border:1px solid #fed7aa;border-radius:8px;padding:0.8rem 1rem;color:#7c2d12;}
    .quality-excellent {background:#f0fdf4;border:1px solid #86efac;border-radius:8px;padding:0.8rem 1rem;color:#14532d;}
    .quality-good {background:#eff6ff;border:1px solid #93c5fd;border-radius:8px;padding:0.8rem 1rem;color:#1e3a5f;}
    .quality-moderate {background:#fffbeb;border:1px solid #fcd34d;border-radius:8px;padding:0.8rem 1rem;color:#78350f;}
    .quality-poor {background:#fef2f2;border:1px solid #fca5a5;border-radius:8px;padding:0.8rem 1rem;color:#7f1d1d;}
    </style>
    """,
    unsafe_allow_html=True,
)


def uploaded_text(uploaded) -> tuple[str | None, str | None]:
    if uploaded is None:
        return None, None
    return uploaded.getvalue().decode("utf-8", errors="replace"), uploaded.name


def slug(text: str) -> str:
    clean = re.sub(r"[^a-zA-Z0-9_]+", "_", text.strip().lower()).strip("_")
    return clean or "custom_component"


def parse_nuclei_text(text: str) -> list[Nucleus]:
    nuclei: list[Nucleus] = []
    if not text or not str(text).strip():
        return nuclei
    for item in re.split(r"[;,]", str(text)):
        item = item.strip()
        if not item:
            continue
        match = re.match(r"([A-Za-z0-9]+)\s*[:=]\s*([0-9.]+)", item)
        if not match:
            continue
        isotope, a_value = match.groups()
        nuclei.append(Nucleus(isotope=isotope, A_mT=float(a_value), label=isotope))
    return nuclei


def custom_components_from_table(table: pd.DataFrame) -> list[SpinComponent]:
    components: list[SpinComponent] = []
    for _, row in table.iterrows():
        name = str(row.get("name", "")).strip()
        if not name:
            continue
        cid = "custom_" + slug(name)
        components.append(
            SpinComponent(
                component_id=cid,
                display_name=name,
                radical_assignment=str(row.get("assignment", name)),
                category=str(row.get("category", "custom")),
                g=float(row.get("g", 2.003)),
                g_bounds=(float(row.get("g_min", 1.95)), float(row.get("g_max", 2.20))),
                nuclei=parse_nuclei_text(str(row.get("nuclei", ""))),
                linewidth_mT=float(row.get("linewidth_mT", 0.2)),
                linewidth_bounds=(float(row.get("lw_min", 0.01)), float(row.get("lw_max", 10.0))),
                eta=float(row.get("eta", 0.5)),
                weight=float(row.get("weight", 0.3)),
                interpretation="User-defined component.",
                warning="Custom component: verify g, hyperfine, and linewidth bounds against literature or controls.",
            )
        )
    return components


def apply_component_edits(selected_ids: list[str], edited: pd.DataFrame, custom_components: list[SpinComponent]) -> list[SpinComponent]:
    library = default_components()
    custom_by_id = {component.component_id: component for component in custom_components}
    components: list[SpinComponent] = []
    for cid in selected_ids:
        if cid in library:
            components.append(library[cid].clone())
        elif cid in custom_by_id:
            components.append(custom_by_id[cid].clone())
    by_id = {component.component_id: component for component in components}
    for _, row in edited.iterrows():
        cid = row.get("component ID")
        if cid in by_id:
            comp = by_id[cid]
            comp.g = float(row.get("g", comp.g))
            comp.linewidth_mT = float(row.get("linewidth mT", comp.linewidth_mT))
            comp.eta = float(row.get("eta", comp.eta))
            comp.weight = float(row.get("weight", comp.weight))
    return components


def anisotropic_editor_frame(components: list[SpinComponent]) -> pd.DataFrame:
    """Build the editable anisotropic-parameter table for the selected components."""
    rows = []
    for comp in components:
        gx, gy, gz = comp.g_principal()
        rows.append({
            "component ID": comp.component_id,
            "S (spin)": float(comp.spin_S),
            "gx": round(gx, 5),
            "gy": round(gy, 5),
            "gz": round(gz, 5),
            "D (MHz)": float(comp.D_MHz),
            "E (MHz)": float(comp.E_MHz),
            "mode": comp.mode,
        })
    cols = ["component ID", "S (spin)", "gx", "gy", "gz", "D (MHz)", "E (MHz)", "mode"]
    return pd.DataFrame(rows, columns=cols)


def apply_anisotropic_edits(components: list[SpinComponent], edited: pd.DataFrame) -> list[SpinComponent]:
    """Apply anisotropic (g-tensor, S, D, E, mode) edits to the selected components."""
    by_id = {c.component_id: c for c in components}
    for _, row in edited.iterrows():
        cid = row.get("component ID")
        comp = by_id.get(cid)
        if comp is None:
            continue
        try:
            comp.spin_S = float(row.get("S (spin)", comp.spin_S))
            gx = float(row.get("gx")); gy = float(row.get("gy")); gz = float(row.get("gz"))
            comp.g_tensor = (gx, gy, gz)
            comp.g = (gx + gy + gz) / 3.0
            comp.D_MHz = float(row.get("D (MHz)", comp.D_MHz))
            comp.E_MHz = float(row.get("E (MHz)", comp.E_MHz))
            mode = str(row.get("mode", comp.mode)).strip().lower()
            comp.mode = mode if mode in ("auto", "isotropic", "powder") else "auto"
        except (TypeError, ValueError):
            continue
    return components


# ── Header banner ────────────────────────────────────────────────────────────
if LOGO_BANNER.exists():
    st.image(str(LOGO_BANNER), use_container_width=True)
else:
    st.title("SimEPR")
st.caption("General high-field cw-EPR simulation, fitting, model comparison, and publication-ready export.")
st.markdown(f"<div class='credit-note'>{SHORT_CREDIT}</div>", unsafe_allow_html=True)
st.markdown(
    "<div class='science-note'>SimEPR uses isotropic high-field cw-EPR line simulations. "
    "It is scientifically useful for screening, mixture fitting, and transparent reporting, "
    "but anisotropic powder spectra, tensor-resolved g/A analysis, saturation behavior, and exchange/correlation effects require specialist EPR treatment.</div>",
    unsafe_allow_html=True,
)

with st.sidebar:
    # Logo in sidebar
    if LOGO.exists():
        st.image(str(LOGO), width=120)
    st.header("Experiment metadata")
    project_title = st.text_input("Project title", value="Untitled EPR study")
    sample_name = st.text_input("Sample name", value="Sample 1")
    solvent = st.text_input("Solvent / matrix", value="water, buffer, frozen glass, solid, or custom")
    catalyst_material = st.text_input("Catalyst/material/combination", value="custom catalyst or material")
    atmosphere = st.text_input("Atmosphere / gas", value="N2, Ar, air, O2, vacuum, or custom")
    condition_text = st.text_input("Condition", value="light/dark, temperature, reaction time, pH, potential, dose, etc.")
    additives = st.text_area("Additives / electrolyte / spin probe", value="", height=70)
    sample_class = st.selectbox(
        "Sample class",
        [
            "unknown/general",
            "organic radical",
            "nitroxide/spin label",
            "ROS/spin trap",
            "transition metal",
            "solid defect/broad signal",
        ],
        index=0,
    )
    pbn_used = st.checkbox("PBN or DMSO-specific spin trapping", value=False)
    st.divider()
    st.header("Data")
    uploaded = st.file_uploader("Upload EPR file", type=["asc", "txt", "dat", "csv"])
    demo_choice = st.selectbox("Or load demo", ["none"] + list(DEMO_CONDITIONS.keys()), index=0)
    mw_freq = st.number_input("Microwave frequency / GHz", min_value=1.0, max_value=300.0, value=DEFAULT_MW_FREQUENCY_GHZ, step=0.01,
                              help="X-band ≈ 9.4 GHz, Q-band ≈ 34 GHz, W-band ≈ 94 GHz. Multifrequency-aware: anisotropic patterns scale with frequency.")
    field_unit = st.selectbox("Field unit", ["Auto", "Gauss", "mT"], index=0)
    advanced = st.checkbox("Advanced fitting controls", value=False)
    st.divider()
    st.header("Powder engine")
    powder_quality = st.select_slider(
        "Orientation accuracy",
        options=["Fast (600)", "Standard (1000)", "High (2000)", "Very high (4000)"],
        value="Standard (1000)",
        help="Orientations for powder averaging of anisotropic / high-spin components. More = smoother, slower.",
    )
    n_orient = {"Fast (600)": 600, "Standard (1000)": 1000, "High (2000)": 2000, "Very high (4000)": 4000}[powder_quality]

file_text, filename = uploaded_text(uploaded)
if file_text is None and demo_choice != "none":
    demo_condition, demo_preset, demo_weights = DEMO_CONDITIONS[demo_choice]
    file_text = generate_demo_text(demo_condition, demo_preset, demo_weights)
    filename = demo_choice

context = ExperimentContext(
    analysis_mode="General EPR fitting",
    sample_class=sample_class,
    condition=condition_text,
    catalyst=catalyst_material,
    solvent=solvent,
    pbn_used=pbn_used,
)
suggested = suggest_models(context)

parsed = None
prep = None
parse_error = None
if file_text:
    try:
        parsed_initial = parse_epr_text(file_text, filename=filename, field_unit=field_unit, mw_frequency_override=mw_freq)
    except Exception as exc:  # noqa: BLE001
        parsed_initial = None
        parse_error = str(exc)
else:
    parsed_initial = None

# Apply any programmatic field-shift written to _pending_field_shift BEFORE the
# slider widget (key="field_shift_mT") is instantiated.  Streamlit forbids writing
# to a widget-bound key after the widget renders, so all button handlers stage their
# update here and call st.rerun() to let this block commit it first.
if "_pending_field_shift" in st.session_state:
    st.session_state["field_shift_mT"] = st.session_state.pop("_pending_field_shift")

tabs = st.tabs(["Import", "Metadata", "Preprocess", "Model builder", "Fit", "Compare", "Export", "White paper / citation"])

with tabs[0]:
    st.subheader("Import")
    with st.expander("📖 Scientific background — cw-EPR fundamentals", expanded=False):
        st.markdown(SCIENCE_IMPORT)
    if parse_error:
        st.error(parse_error)
    if parsed_initial is None:
        st.info("Upload an EPR text/ASC/CSV file or load a demo.")
    else:
        col_opts = list(range(parsed_initial.detected_columns))
        c1, c2, c3, c4 = st.columns(4)
        field_col = c1.selectbox("Field column", col_opts, index=min(parsed_initial.field_col, len(col_opts) - 1))
        intensity_col = c2.selectbox("Intensity column", col_opts, index=min(1, len(col_opts) - 1))
        manual_unit = c3.selectbox("Manual field unit", ["Auto", "Gauss", "mT"], index=["Auto", "Gauss", "mT"].index(field_unit))
        manual_mw = c4.number_input("Manual microwave GHz", min_value=1.0, max_value=300.0, value=float(mw_freq), step=0.01)
        parsed = parse_epr_text(
            file_text,
            filename=filename,
            field_col=int(field_col),
            intensity_col=int(intensity_col),
            field_unit=manual_unit,
            mw_frequency_override=float(manual_mw),
        )
        context.metadata = parsed.metadata
        suggested = suggest_models(context)
        m1, m2, m3, m4, m5 = st.columns(5)
        m1.metric("File", parsed.filename or "uploaded")
        m2.metric("Rows", f"{parsed.numeric_rows:,}")
        m3.metric("Columns", parsed.detected_columns)
        m4.metric("Field unit", parsed.detected_field_unit)
        m5.metric("MW freq", f"{parsed.metadata.get('microwave_frequency_GHz', manual_mw):.4g} GHz")
        for warning in parsed.warnings:
            st.warning(warning)
        st.dataframe(pd.DataFrame(metadata_table(parsed.metadata)), width="stretch")
        with st.expander("Header text"):
            st.text(parsed.header_text or "(no non-numeric header detected)")
        st.dataframe(parsed.dataframe.head(500), width="stretch")
        st.plotly_chart(spectrum_figure(parsed.dataframe["Field_mT"], {"raw": parsed.dataframe["Intensity_raw"]}, "Raw EPR spectrum"), width="stretch")

with tabs[1]:
    st.subheader("Metadata and suggested models")
    meta_rows = [
        ("Project", project_title),
        ("Sample", sample_name),
        ("Solvent/matrix", solvent),
        ("Catalyst/material", catalyst_material),
        ("Atmosphere/gas", atmosphere),
        ("Condition", condition_text),
        ("Additives/spin probe", additives),
        ("Sample class", sample_class),
    ]
    st.dataframe(pd.DataFrame(meta_rows, columns=["field", "value"]), width="stretch")
    if pbn_used and "dmso" in (solvent + additives).lower():
        st.markdown(
            "<div class='warning-note'>PBN/DMSO was indicated. DMSO-derived carbon radical adducts must be considered before assigning other spin-trapped species.</div>",
            unsafe_allow_html=True,
        )
    cards = []
    for model in suggested["suggested_models"] + [m for m in MODEL_PRESETS if m not in suggested["suggested_models"]]:
        cards.append(
            {
                "model": model,
                "suggested": model in suggested["suggested_models"],
                "description": MODEL_DESCRIPTIONS.get(model, ""),
                "components": ", ".join(MODEL_PRESETS[model]),
            }
        )
    st.dataframe(pd.DataFrame(cards), width="stretch")

with tabs[2]:
    st.subheader("Preprocess")
    with st.expander("📖 Scientific background — baseline correction, normalisation & field alignment", expanded=False):
        st.markdown(SCIENCE_PREPROCESS)
    if parsed is None:
        st.info("Import a spectrum first.")
    else:
        raw_field = parsed.dataframe["Field_mT"].to_numpy()
        raw_intensity = parsed.dataframe["Intensity_raw"].to_numpy()

        # Reset shift when a new file is loaded
        _file_key = filename or ""
        if st.session_state.get("_preproc_file") != _file_key:
            st.session_state["field_shift_mT"] = 0.0
            st.session_state["_preproc_file"] = _file_key
            st.session_state.pop("sim_preview", None)
        if "field_shift_mT" not in st.session_state:
            st.session_state["field_shift_mT"] = 0.0

        field_shift = float(st.session_state["field_shift_mT"])
        shifted_min = float(raw_field.min()) + field_shift
        shifted_max = float(raw_field.max()) + field_shift

        pc1, pc2, pc3, pc4 = st.columns(4)
        crop_min, crop_max = pc1.slider("Fitting window / mT", shifted_min, shifted_max, (shifted_min, shifted_max))
        baseline_method = pc2.selectbox("Baseline correction", ["none", "constant", "linear edge", "polynomial edge"], index=2)
        poly_order = pc2.number_input("Polynomial order", min_value=1, max_value=5, value=2)
        norm_method = pc3.selectbox("Normalization", ["none", "max absolute", "peak-to-peak", "area"], index=1)
        invert = pc3.checkbox("Invert intensity", value=False)
        smooth = pc4.checkbox("Display smoothing", value=False)

        prep = preprocess_spectrum(
            raw_field,
            raw_intensity,
            crop_min=crop_min,
            crop_max=crop_max,
            field_shift_mT=field_shift,
            baseline_method=baseline_method,
            polynomial_order=int(poly_order),
            normalization=norm_method,
            smooth_display=smooth,
            invert=invert,
        )

        preprocess_title = "Preprocessing" + (f"  ·  field shift {field_shift:+.3f} mT" if field_shift else "")
        st.plotly_chart(
            spectrum_figure(
                prep.field_mT,
                {"raw window": prep.raw, "baseline corrected": prep.corrected, "processed": prep.processed, "display": prep.display},
                preprocess_title,
            ),
            width="stretch",
        )

        # ── Field alignment ──────────────────────────────────────────────────
        with st.expander("Field alignment — shift experimental axis to match model peak", expanded=(field_shift != 0.0)):
            st.markdown(
                "Compensate for spectrometer field-axis calibration offsets. "
                "Build a model in the **Model builder** tab, then click **Auto-align** to "
                "compute the optimal shift via cross-correlation. "
                "Fine-tune manually with the slider."
            )
            acol1, acol2, acol3 = st.columns([4, 2, 1])
            acol1.slider(
                "Field shift / mT",
                min_value=-10.0,
                max_value=10.0,
                step=0.01,
                format="%.2f",
                key="field_shift_mT",
                help="Shifts the experimental field axis. Positive = axis moves up (spectrum shifts left on plot).",
            )
            sim_preview = st.session_state.get("sim_preview")
            auto_disabled = sim_preview is None
            auto_help = (
                "Compute the optimal field shift by cross-correlating the preprocessed experimental "
                "spectrum with the current model simulation."
                if not auto_disabled
                else "Build a model in the Model builder tab first."
            )
            if acol2.button("Auto-align with model", disabled=auto_disabled, help=auto_help):
                sim_field_prev, sim_total_prev = sim_preview  # type: ignore[misc]
                if len(prep.field_mT) > 1 and len(sim_total_prev) > 1:
                    from epr_simfit.preprocessing import auto_align_spectra  # lazy import avoids hot-reload cache issues
                    sim_on_grid = np.interp(prep.field_mT, sim_field_prev, sim_total_prev)
                    delta = auto_align_spectra(prep.processed, sim_on_grid, prep.field_mT)
                    st.session_state["_pending_field_shift"] = round(field_shift + delta, 2)
                    st.rerun()
            if acol3.button("Reset", help="Reset field shift to zero"):
                st.session_state["_pending_field_shift"] = 0.0
                st.rerun()
            if field_shift:
                st.info(f"Active shift: **{field_shift:+.3f} mT** — the experimental field axis has been offset by this amount.")
            elif auto_disabled:
                st.caption("Build a model in the Model builder tab, then click Auto-align to set the shift automatically.")

            # ── Live alignment preview ──────────────────────────────────────
            _sp = st.session_state.get("sim_preview")
            if _sp is not None:
                _sf, _st2 = _sp
                _sim_live = np.interp(prep.field_mT, _sf, _st2)
                st.plotly_chart(
                    spectrum_figure(
                        prep.field_mT,
                        {"experimental (shifted)": prep.processed, "simulation": _sim_live},
                        f"Live alignment preview  ·  shift = {field_shift:+.3f} mT",
                    ),
                    use_container_width=True,
                )

with tabs[3]:
    st.subheader("Model builder")
    with st.expander("📖 Scientific background — isotropic cw-EPR simulation", expanded=False):
        st.markdown(SCIENCE_MODELBUILDER)
    st.dataframe(pd.DataFrame(component_table()), width="stretch")
    preset_names = list(MODEL_PRESETS.keys()) + ["custom"]
    default_preset = suggested.get("recommended_model", "G0_single_line")
    preset = st.selectbox("Preset", preset_names, index=preset_names.index(default_preset) if default_preset in preset_names else 0)
    default_ids = MODEL_PRESETS.get(preset, suggested["suggested_models"][:1])
    st.caption("Custom nuclei syntax: `14N:1.55; 1H:0.30; 63Cu:8.5` where hyperfine values are in mT.")
    custom_default = pd.DataFrame(
        [
            {
                "name": "",
                "assignment": "",
                "category": "custom",
                "g": 2.003,
                "g_min": 1.95,
                "g_max": 2.20,
                "nuclei": "",
                "linewidth_mT": 0.20,
                "lw_min": 0.01,
                "lw_max": 10.0,
                "eta": 0.5,
                "weight": 0.3,
            }
        ]
    )
    custom_table = st.data_editor(custom_default, num_rows="dynamic", width="stretch")
    custom_components = custom_components_from_table(custom_table)
    component_ids = list(default_components().keys()) + [component.component_id for component in custom_components]
    selected_ids = st.multiselect("Components to simulate/fit", component_ids, default=[cid for cid in default_ids if cid in component_ids])
    rows = []
    all_by_id = {**default_components(), **{component.component_id: component for component in custom_components}}
    for cid in selected_ids:
        rows.append(all_by_id[cid].to_table_row())
    edit_cols = ["component ID", "name", "category", "g", "linewidth mT", "eta", "weight", "interpretation", "warning"]
    edited = st.data_editor(pd.DataFrame(rows, columns=edit_cols), disabled=["component ID", "name", "category", "interpretation", "warning"], width="stretch")
    selected_components = apply_component_edits(selected_ids, edited, custom_components)

    # ── Anisotropic / high-spin parameters (tensor g, S, zero-field splitting) ──
    with st.expander("⚛ Anisotropic / high-spin parameters (g-tensor, spin S, zero-field splitting)", expanded=False):
        st.markdown(
            "Set the **g-tensor** principal values (gx, gy, gz), electron **spin S**, and "
            "**zero-field splitting** D, E (MHz) for rigid-limit, powder, frozen-solution, and "
            "high-spin systems. When gx≠gy≠gz, S>½, or D/E≠0, the component is solved by full "
            "spin-Hamiltonian diagonalisation with powder averaging. "
            "`mode`: *auto* (detect), *isotropic* (force fast path), *powder* (force averaging)."
        )
        if selected_components:
            aniso_frame = anisotropic_editor_frame(selected_components)
            aniso_edited = st.data_editor(
                aniso_frame, disabled=["component ID"], width="stretch", key="aniso_editor",
                column_config={
                    "S (spin)": st.column_config.NumberColumn(help="Electron spin: 0.5, 1, 1.5, 2, 2.5 ...", min_value=0.5, max_value=3.5, step=0.5),
                    "D (MHz)": st.column_config.NumberColumn(help="Axial zero-field splitting (S>1/2 only)"),
                    "E (MHz)": st.column_config.NumberColumn(help="Rhombic zero-field splitting (|E/D| ≤ 1/3)"),
                },
            )
            selected_components = apply_anisotropic_edits(selected_components, aniso_edited)
            _any_aniso = any(c.is_anisotropic() for c in selected_components)
            if _any_aniso:
                st.info(f"Powder averaging active for anisotropic/high-spin components ({n_orient} orientations). "
                        "Increase orientation accuracy in the sidebar for smoother patterns.")
        else:
            st.caption("Select components above to edit their anisotropic parameters.")

    if parsed is not None and prep is not None and selected_components:
        weights = {component.component_id: component.weight for component in selected_components}
        _need_powder = any(c.is_anisotropic() for c in selected_components)
        if _need_powder:
            with st.spinner(f"Simulating powder pattern ({n_orient} orientations)..."):
                sim_total, sim_curves = simulate_model(prep.field_mT, selected_components, weights=weights, mw_frequency_GHz=mw_freq, n_orientations=n_orient)
        else:
            sim_total, sim_curves = simulate_model(prep.field_mT, selected_components, weights=weights, mw_frequency_GHz=mw_freq, n_orientations=n_orient)
        # Cache for auto-alignment in the Preprocess tab
        st.session_state["sim_preview"] = (prep.field_mT.copy(), sim_total.copy())
        st.session_state["selected_components"] = selected_components
        st.session_state["n_orient"] = n_orient
        traces = {"experimental processed": prep.processed, "simulation total": sim_total}
        traces.update({cid: weights.get(cid, 1.0) * curve for cid, curve in sim_curves.items()})
        st.plotly_chart(spectrum_figure(prep.field_mT, traces, "Manual simulation"), width="stretch")

with tabs[4]:
    st.subheader("Fit")
    with st.expander("📖 Scientific background — fitting algorithm, parameters & statistics", expanded=False):
        st.markdown(SCIENCE_FIT)
    if parsed is None or prep is None or not selected_components:
        st.info("Import, preprocess, and select components first.")
    else:
        fc1, fc2, fc3 = st.columns(3)
        mode = fc1.selectbox("Fit mode", ["weights only", "weights + linewidths", "weights + linewidths + g"], index=0)
        baseline_order = fc2.selectbox("Baseline in fit", [0, 1], index=0)
        max_eval = fc3.number_input("Max evaluations", min_value=100, max_value=5000, value=700, step=100)

        # ── Spectral masking ──────────────────────────────────────────────────
        with st.expander("Spectral masking — exclude field regions from fit optimization", expanded=False):
            st.markdown(
                "Define field ranges to **exclude** from the least-squares optimization. "
                "Useful when peaks in those regions are artefacts, solvent lines, or belong to an unmodelled species. "
                "Excluded regions are still shown in the overlay plot so you can inspect them."
            )
            _mask_default = pd.DataFrame({"From (mT)": pd.Series([], dtype=float), "To (mT)": pd.Series([], dtype=float)})
            _mask_state = st.session_state.get("fit_masks_df", _mask_default.to_dict("records"))
            _mask_edited = st.data_editor(
                pd.DataFrame(_mask_state) if _mask_state else _mask_default,
                num_rows="dynamic",
                use_container_width=True,
                key="mask_editor",
                column_config={
                    "From (mT)": st.column_config.NumberColumn("From (mT)", format="%.2f"),
                    "To (mT)": st.column_config.NumberColumn("To (mT)", format="%.2f"),
                },
            )
            st.session_state["fit_masks_df"] = _mask_edited.to_dict("records")
            _fit_masks = [
                (float(r["From (mT)"]), float(r["To (mT)"]))
                for r in st.session_state["fit_masks_df"]
                if r.get("From (mT)") is not None and r.get("To (mT)") is not None
                and float(r["From (mT)"]) < float(r["To (mT)"])
            ]

        # Build boolean keep-mask for fitting
        _fit_keep = np.ones(len(prep.field_mT), dtype=bool)
        for _mlo, _mhi in _fit_masks:
            _fit_keep &= ~((prep.field_mT >= _mlo) & (prep.field_mT <= _mhi))
        _n_excl = int((~_fit_keep).sum())
        if _n_excl:
            st.caption(f"ℹ {_n_excl} of {len(prep.field_mT)} data points are excluded from the fit by spectral masks.")

        # ── Run fit ───────────────────────────────────────────────────────────
        _auto_refit = st.session_state.pop("_auto_refit", False)
        if st.button("Run fit", type="primary") or _auto_refit:
            _use_field = prep.field_mT[_fit_keep]
            _use_exp   = prep.processed[_fit_keep]
            _spinner_msg = "Re-fitting with aligned field..." if _auto_refit else (
                f"Fitting on {len(_use_field)} points ({_n_excl} masked)..." if _n_excl else "Fitting..."
            )
            with st.spinner(_spinner_msg):
                fit = fit_spectrum(
                    _use_field,
                    _use_exp,
                    preset=preset,
                    components=selected_components,
                    mw_frequency_GHz=mw_freq,
                    mode=mode,
                    baseline_order=int(baseline_order),
                    max_nfev=int(max_eval),
                    n_orientations=st.session_state.get("n_orient", 1000),
                )
                # Re-evaluate on FULL grid for display (masked fit gives parameters only)
                if _n_excl > 0:
                    from dataclasses import replace as _dc_replace
                    _full_sim, _full_curves = simulate_model(
                        prep.field_mT, fit.components, fit.weights, mw_freq, fit.baseline0, fit.baseline1
                    )
                    fit = _dc_replace(
                        fit,
                        field_mT=prep.field_mT,
                        experimental=prep.processed,
                        fit_total=_full_sim,
                        residual=prep.processed - _full_sim,
                        component_curves=_full_curves,
                    )
            st.session_state["fit"] = fit

        fit = st.session_state.get("fit")
        if fit:
            # ── Fit overlay ──────────────────────────────────────────────────
            _overlay_traces = {"experimental": fit.experimental, "fit": fit.fit_total}
            if _fit_masks:
                _mask_show = np.full_like(fit.experimental, np.nan)
                for _mlo, _mhi in _fit_masks:
                    _mi = (fit.field_mT >= _mlo) & (fit.field_mT <= _mhi)
                    _mask_show[_mi] = fit.experimental[_mi]
                _overlay_traces["excluded (masked)"] = _mask_show
            _overlay_fig = spectrum_figure(fit.field_mT, _overlay_traces, "Fit overlay")
            st.plotly_chart(_overlay_fig, use_container_width=True)

            # ── Auto-align from fit overlay ──────────────────────────────────
            _cur_shift = float(st.session_state.get("field_shift_mT", 0.0))
            _oa1, _oa2 = st.columns([2, 5])
            if _oa1.button("Auto-align peak to fit",
                           help="Shift experimental axis to the fitted peak, then refit."):
                from epr_simfit.preprocessing import auto_align_spectra
                _delta = auto_align_spectra(fit.experimental, fit.fit_total, fit.field_mT)
                if abs(_delta) > 0.005:
                    st.session_state["_pending_field_shift"] = round(_cur_shift + _delta, 2)
                    st.session_state["_auto_refit"] = True
                    st.rerun()
                else:
                    st.toast("Peaks are already aligned.", icon="✅")
            _oa2.caption(f"Active field shift: **{_cur_shift:+.2f} mT** — adjust in Preprocess → Field alignment.")

            # ── Residual plot + noise-floor band ─────────────────────────────
            st.markdown("**Residual  (experimental − fit)**")
            _rc1, _rc2 = st.columns([4, 1])
            _res_pct = _rc2.slider("Noise floor %", 0, 50, 0, 1,
                help="Gray band at ±N% of peak residual. Features outside the band are likely real unmodelled signal.")
            _res_fig = residual_figure(fit.field_mT, fit.residual)
            if _res_pct > 0:
                _abs_max = float(np.nanmax(np.abs(fit.residual))) if fit.residual.size else 1.0
                _thresh  = _abs_max * _res_pct / 100.0
                _res_fig.add_hrect(y0=-_thresh, y1=_thresh,
                                   fillcolor="rgba(180,180,180,0.25)",
                                   line_width=1, line_color="rgba(150,150,150,0.5)",
                                   annotation_text=f"±{_res_pct}% noise floor",
                                   annotation_position="top right")
            _rc1.plotly_chart(_res_fig, use_container_width=True)

            # ── Goodness of fit ───────────────────────────────────────────────
            st.subheader("Goodness of fit")
            fq = assess_fit_quality(fit.metrics)
            qc = st.columns(6)
            qc[0].metric("Status",       f"{fq.icon} {fq.label}")
            qc[1].metric("R²",           f"{fit.metrics['R2']:.4f}")
            qc[2].metric("Norm. RMSE",   f"{fit.metrics['normalized RMSE']:.4f}")
            qc[3].metric("RMSE",         f"{fit.metrics['RMSE']:.4g}")
            qc[4].metric("AIC",          f"{fit.metrics['AIC']:.1f}")
            qc[5].metric("BIC",          f"{fit.metrics['BIC']:.1f}")
            st.markdown(f"**{fq.r2_note}**")
            st.caption(fq.nrmse_note)
            st.info(fq.overall)

            # ── Improvement suggestions ───────────────────────────────────────
            _suggestions = suggest_fit_improvements(fit, fit_mode=mode)
            if _suggestions:
                with st.expander(
                    f"💡 Suggestions to improve this fit  ({len(_suggestions)} item{'s' if len(_suggestions)!=1 else ''})",
                    expanded=any(s.priority == "High" for s in _suggestions),
                ):
                    for _sg in _suggestions:
                        _col_icon = {"High": "🔴", "Medium": "🟡", "Low": "🔵", "Info": "✅"}.get(_sg.priority, "ℹ️")
                        st.markdown(f"**{_col_icon} [{_sg.priority}] {_sg.title}**")
                        st.caption(f"*Why:* {_sg.reason}")
                        st.caption(f"*What to do:* {_sg.action}")
                        st.divider()

            # ── Component decomposition ───────────────────────────────────────
            _comp_fig = component_figure(fit.field_mT, fit.component_curves, fit.weights)
            st.plotly_chart(_comp_fig, use_container_width=True)

            # ── Extend model and refit ────────────────────────────────────────
            _fitted_ids = [c.component_id for c in fit.components]
            _extendable = [cid for cid in component_ids if cid not in _fitted_ids]
            with st.expander("➕ Extend model — add components and refit", expanded=False):
                st.caption(
                    "Pick additional components to add to the fitted model and immediately refit. "
                    "Useful when the residual shows unaccounted peaks."
                )
                _extra_ids = st.multiselect("Components to add", _extendable, key="extend_component_ids")
                if _extra_ids:
                    _extra_comps = [all_by_id[cid].clone() for cid in _extra_ids if cid in all_by_id]
                    _extended_comps = fit.components + _extra_comps
                    _ex1, _ex2 = st.columns([2, 5])
                    if _ex1.button("Refit with extended model", type="primary", key="btn_refit_extended"):
                        with st.spinner(f"Refitting with {len(_extended_comps)} components..."):
                            _ext_fit = fit_spectrum(
                                prep.field_mT[_fit_keep], prep.processed[_fit_keep],
                                preset=preset, components=_extended_comps,
                                mw_frequency_GHz=mw_freq, mode=mode,
                                baseline_order=int(baseline_order), max_nfev=int(max_eval),
                                n_orientations=st.session_state.get("n_orient", 1000),
                            )
                            if _n_excl > 0:
                                from dataclasses import replace as _dc_replace
                                _fs2, _fc2 = simulate_model(
                                    prep.field_mT, _ext_fit.components, _ext_fit.weights,
                                    mw_freq, _ext_fit.baseline0, _ext_fit.baseline1,
                                )
                                _ext_fit = _dc_replace(_ext_fit,
                                    field_mT=prep.field_mT, experimental=prep.processed,
                                    fit_total=_fs2, residual=prep.processed - _fs2,
                                    component_curves=_fc2)
                        st.session_state["fit"] = _ext_fit
                        st.rerun()
                    _ex2.caption(
                        f"Will add: {', '.join(_extra_ids)}. "
                        f"New model has {len(_extended_comps)} components vs current {len(fit.components)}."
                    )

            # ── Save to comparison ────────────────────────────────────────────
            _fs_val     = float(st.session_state.get("field_shift_mT", 0.0))
            _save_label = f"{fit.preset}  shift={_fs_val:+.2f}mT  R²={fit.metrics['R2']:.3f}"
            sc1, sc2 = st.columns([2, 5])
            if sc1.button("💾 Save fit to comparison", help="Add to the Compare tab for side-by-side model selection."):
                if "saved_fits" not in st.session_state:
                    st.session_state["saved_fits"] = {}
                _n = len(st.session_state["saved_fits"]) + 1
                st.session_state["saved_fits"][f"[{_n}] {_save_label}"] = fit
                st.toast(f"Saved as '[{_n}] {_save_label}'", icon="💾")
            sc2.caption(f"Will be saved as: **{_save_label}**")

            # ── Export plots ──────────────────────────────────────────────────
            with st.expander("📥 Export plots", expanded=False):
                st.caption("Download interactive HTML files — open in any browser, zoom, pan, hover for values.")
                _ep1, _ep2, _ep3 = st.columns(3)
                _ep1.download_button(
                    "Fit overlay (HTML)",
                    data=_overlay_fig.to_html(include_plotlyjs="cdn").encode(),
                    file_name="SimEPR_fit_overlay.html", mime="text/html",
                )
                _ep2.download_button(
                    "Residual (HTML)",
                    data=_res_fig.to_html(include_plotlyjs="cdn").encode(),
                    file_name="SimEPR_residual.html", mime="text/html",
                )
                _ep3.download_button(
                    "Components (HTML)",
                    data=_comp_fig.to_html(include_plotlyjs="cdn").encode(),
                    file_name="SimEPR_components.html", mime="text/html",
                )

            # ── Detected intermediates ────────────────────────────────────────
            st.subheader("Detected paramagnetic intermediates / species")
            _it1, _it2 = st.columns([3, 1])
            _interm_thresh = _it2.slider("Min fraction %", 1, 30, 3, 1,
                help="Only show components contributing at least this fraction of total spectral intensity.")
            _it1.caption(
                f"Components ≥ {_interm_thresh}% shown. "
                "Assignments are candidate identifications — validate with controls and isotope labelling."
            )
            intermediates = suggest_intermediates(fit, threshold_pct=float(_interm_thresh))
            if intermediates:
                st.dataframe(intermediates_dataframe(intermediates), use_container_width=True)
                for d in intermediates:
                    with st.expander(f"{d.rank}. {d.name}  —  {d.fraction_pct:.1f}%  ({d.confidence})", expanded=False):
                        st.markdown(f"**Assignment:** {d.assignment}  |  **Category:** {d.category}")
                        st.markdown(f"**g-value:** `{d.g:.5f}`  |  **ΔBpp:** `{d.linewidth_mT:.4f} mT`  |  **Hyperfine:** {d.nuclei_str}")
                        if d.interpretation:
                            st.info(d.interpretation)
                        if d.warning:
                            st.warning(d.warning)
            else:
                st.info(f"No components exceed {_interm_thresh}%. Lower the slider or add more components.")

            # ── Publication-ready parameters ──────────────────────────────────
            st.subheader("Publication-ready parameters")
            pub_df = publication_parameters_table(fit)
            st.dataframe(pub_df, use_container_width=True)
            with st.expander("Detailed fitted parameter table (raw)", expanded=False):
                st.dataframe(fit.parameters, use_container_width=True)
                st.dataframe(fit.component_fractions, use_container_width=True)

            # ── Methods paragraph ─────────────────────────────────────────────
            st.subheader("Publication methods paragraph")
            st.caption("Copy into your manuscript methods or supplementary information.")
            _methods_text = publication_methods_paragraph(fit, mw_freq, float(st.session_state.get("field_shift_mT", 0.0)))
            st.text_area("Methods paragraph", _methods_text, height=220, label_visibility="collapsed")
            st.download_button(
                "⬇ Download publication parameters (CSV)",
                data=pub_df.to_csv(index=False).encode(),
                file_name="SimEPR_publication_parameters.csv", mime="text/csv",
            )
            st.warning("Fit quality is not chemical proof. Validate with standards, controls, isotope/substitution tests, and chemistry-specific constraints.")

with tabs[5]:
    st.subheader("Compare")
    st.caption(
        "Save fits using **💾 Save fit to comparison** in the Fit tab after each run "
        "(different models, modes, alignments). This tab compares them side-by-side."
    )

    _saved_fits: dict = st.session_state.get("saved_fits", {})

    if not _saved_fits:
        st.info(
            "No fits saved yet.  \n"
            "1. Go to the **Fit** tab and run a fit.  \n"
            "2. Click **💾 Save fit to comparison**.  \n"
            "3. Repeat with different models or settings.  \n"
            "4. Return here to compare them by BIC / AIC / R²."
        )
    else:
        # ── Metrics table ─────────────────────────────────────────────────────
        _sf_rows = []
        for _lbl, _sf in _saved_fits.items():
            _sf_rows.append({"Model": _lbl, "n params": _sf.n_parameters, **{k: v for k, v in _sf.metrics.items()}})
        _sf_df = pd.DataFrame(_sf_rows)
        _sf_df["ΔBIC"] = (_sf_df["BIC"] - _sf_df["BIC"].min()).round(2)
        _sf_df["ΔAIC"] = (_sf_df["AIC"] - _sf_df["AIC"].min()).round(2)
        _sf_df = _sf_df.sort_values("BIC").reset_index(drop=True)

        # Highlight best row
        _best_label = _sf_df.iloc[0]["Model"]
        st.markdown(f"**Best model by BIC: {_best_label}** — ΔBIC > 10 is strong evidence for preference.")
        st.dataframe(_sf_df, use_container_width=True)

        # ── BIC bar chart ─────────────────────────────────────────────────────
        _cmp_for_bar = _sf_df.rename(columns={"Model": "model"})[["model", "BIC", "AIC", "R2", "ΔBIC", "ΔAIC"]]
        _metric_sel = st.radio("Plot metric", ["BIC", "AIC", "R2"], horizontal=True, index=0)
        st.plotly_chart(comparison_bar_figure(_cmp_for_bar, _metric_sel), use_container_width=True)

        # ── Overlay all saved fits on experimental spectrum ───────────────────
        if parsed is not None and prep is not None:
            with st.expander("Overlay all saved fits on experimental spectrum", expanded=True):
                _ov_traces = {"experimental": prep.processed}
                for _lbl, _sf in _saved_fits.items():
                    _ov_traces[_lbl] = np.interp(prep.field_mT, _sf.field_mT, _sf.fit_total)
                _ov_fig = spectrum_figure(prep.field_mT, _ov_traces, "All saved fits vs experimental")
                st.plotly_chart(_ov_fig, use_container_width=True)
                st.download_button(
                    "⬇ Export comparison overlay (HTML)",
                    data=_ov_fig.to_html(include_plotlyjs="cdn").encode(),
                    file_name="SimEPR_comparison_overlay.html", mime="text/html",
                )

        # ── BIC interpretation guide ──────────────────────────────────────────
        with st.expander("ΔBIC interpretation (Kass & Raftery, 1995)", expanded=False):
            st.markdown("""
| ΔBIC | Evidence against higher-BIC model |
|------|----------------------------------|
| 0 – 2 | Not worth more than a bare mention |
| 2 – 6 | Positive (weak) evidence |
| 6 – 10 | Strong evidence |
| > 10 | Very strong evidence |

Lower BIC = statistically preferred. BIC penalises additional parameters more heavily than AIC.
Use model comparison as a guide; chemical validation always takes precedence.
""")

        # ── Download comparison table ──────────────────────────────────────────
        _dc1, _dc2 = st.columns([2, 5])
        _dc1.download_button(
            "⬇ Comparison table (CSV)",
            data=_sf_df.to_csv(index=False).encode(),
            file_name="SimEPR_model_comparison.csv", mime="text/csv",
        )
        if _dc2.button("🗑 Clear all saved fits"):
            st.session_state["saved_fits"] = {}
            st.rerun()

    # ── Batch preset comparison (secondary) ──────────────────────────────────
    with st.expander("Run batch preset comparison on current aligned spectrum (advanced)", expanded=False):
        st.caption(
            "Fits the selected library presets on the current preprocessed + aligned data. "
            "Results appear only here and are NOT added to the saved-fit comparison above."
        )
        if parsed is None or prep is None:
            st.info("Import and preprocess first.")
        else:
            defaults2 = ["G0_single_line"] + [m for m in suggested["suggested_models"] if m != "G0_single_line"][:3]
            compare_presets = st.multiselect("Presets", list(MODEL_PRESETS.keys()), default=list(dict.fromkeys(defaults2)), key="cmp_presets2")
            compare_mode = st.selectbox("Fit mode", ["weights only", "weights + linewidths"], index=0, key="cmp_mode2")
            if st.button("Run batch comparison"):
                with st.spinner("Running batch comparison on aligned spectrum..."):
                    comparison, _ = compare_models(prep.field_mT, prep.processed, compare_presets, mw_frequency_GHz=mw_freq, mode=compare_mode)
                st.session_state["batch_comparison"] = comparison
            _bc = st.session_state.get("batch_comparison")
            if _bc is not None and not _bc.empty:
                st.dataframe(_bc, use_container_width=True)
                st.plotly_chart(comparison_bar_figure(_bc, "BIC"), use_container_width=True)

with tabs[6]:
    st.subheader("Export")
    fit = st.session_state.get("fit")
    comparison = st.session_state.get("comparison")
    processed = None
    if parsed is not None and prep is not None:
        processed = pd.DataFrame({"Field_mT": prep.field_mT, "Intensity_raw": prep.raw, "Intensity_processed": prep.processed})
    evidence_summary = {
        "evidence_class": "GENERAL_EPR_FIT",
        "explanation": "SimEPR reports candidate spectral decompositions using high-field isotropic cw-EPR models.",
        "warnings": ["Validate chemical assignments independently; anisotropic tensor cases require specialist EPR analysis."],
        "recommended_manuscript_wording": "The cw-EPR spectrum was fitted using SimEPR high-field isotropic model components, and candidate assignments were evaluated with model comparison and chemical controls.",
    }
    export_context = {
        "project_title": project_title,
        "sample_name": sample_name,
        "solvent_matrix": solvent,
        "catalyst_material_combination": catalyst_material,
        "atmosphere_gas": atmosphere,
        "condition": condition_text,
        "additives_spin_probe": additives,
        "sample_class": sample_class,
        "software": "SimEPR",
    }
    report_txt = generate_report_text(
        parsed.metadata if parsed else {},
        context,
        prep.settings if prep else {},
        fit.preset if fit else preset if "preset" in locals() else "not fitted",
        fit.parameters if fit else None,
        comparison,
        evidence_summary,
    )
    report_txt = "User-entered experiment metadata\n" + str(export_context) + "\n\n" + report_txt
    st.text_area("Report preview", report_txt, height=320)
    if fit:
        st.dataframe(fit_components_dataframe(fit).head(300), width="stretch")
    if st.button("Prepare SimEPR export package"):
        with st.spinner("Preparing CSVs, fitted metrics, plot data, figures, report, and citation files..."):
            zip_bytes = build_export_zip(
                fit=fit,
                processed=processed,
                comparison=comparison,
                evidence_summary=evidence_summary,
                report_txt=report_txt,
                report_html=generate_report_html(report_txt),
                config=export_context,
                mw_frequency_GHz=mw_freq,
            )
            if WHITE_PAPER.exists() or CITATION.exists():
                from io import BytesIO
                from zipfile import ZIP_DEFLATED, ZipFile

                buf = BytesIO()
                with ZipFile(BytesIO(zip_bytes)) as zin, ZipFile(buf, "w", ZIP_DEFLATED) as zout:
                    for item in zin.infolist():
                        zout.writestr(item, zin.read(item.filename))
                    if WHITE_PAPER.exists():
                        zout.writestr("citation/WHITE_PAPER.md", WHITE_PAPER.read_text(encoding="utf-8"))
                    if CITATION.exists():
                        zout.writestr("citation/CITATION.cff", CITATION.read_text(encoding="utf-8"))
                zip_bytes = buf.getvalue()
            st.session_state["export_zip"] = zip_bytes
    st.caption("Export includes fit metrics, model metrics, all plotted datasets as CSV, figures, reports, citation files, EasySpin script, and ORCA templates where a fit is available.")
    st.download_button("Download SimEPR export ZIP", data=st.session_state.get("export_zip", b""), file_name="SimEPR_export.zip", mime="application/zip", disabled="export_zip" not in st.session_state)

with tabs[7]:
    st.subheader("White paper / citation")
    st.write(ABOUT_TEXT)
    st.markdown("**Developer**")
    st.write(DETAILED_INFO)
    if WHITE_PAPER.exists():
        whitepaper = WHITE_PAPER.read_text(encoding="utf-8")
        st.download_button("Download SimEPR white paper", data=whitepaper, file_name="SimEPR_WHITE_PAPER.md", mime="text/markdown")
        st.markdown(whitepaper)
    if CITATION.exists():
        citation = CITATION.read_text(encoding="utf-8")
        st.download_button("Download CITATION.cff", data=citation, file_name="CITATION.cff", mime="text/plain")
        st.code(citation, language="yaml")

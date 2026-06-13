import pandas as pd

from epr_simfit.model_library import MODEL_PRESETS, default_components
from epr_simfit.model_suggester import ExperimentContext, suggest_models
from epr_simfit.spin_models import Nucleus, generate_hyperfine_lines
from epr_simfit.user_models import components_from_csv, components_from_json, components_to_json


def test_public_general_transition_metal_suggestion():
    out = suggest_models(ExperimentContext(analysis_mode="General EPR fitting", sample_class="transition metal", catalyst="user catalyst"))
    assert out["recommended_model"] == "G4_transition_metal_screen"


def test_supported_public_nuclei_include_mn_and_cu():
    mn_lines, _ = generate_hyperfine_lines([Nucleus(isotope="55Mn", A_mT=9.4)])
    cu_lines, _ = generate_hyperfine_lines([Nucleus(isotope="63Cu", A_mT=8.5)])
    assert len(mn_lines) == 6
    assert len(cu_lines) == 4


def test_co2rr_and_n2rr_routes_have_real_presets():
    co2 = suggest_models(ExperimentContext(analysis_mode="General EPR fitting", sample_class="CO2RR / carbon dioxide reduction"))
    n2 = suggest_models(ExperimentContext(analysis_mode="General EPR fitting", sample_class="N2RR / nitrogen reduction"))
    assert "C1_co2rr_screen" in co2["suggested_models"]
    assert "N1_n2rr_dmso_corrected" in n2["suggested_models"]
    for name in co2["suggested_models"] + n2["suggested_models"]:
        assert name in MODEL_PRESETS


def test_electrocatalysis_route_has_chemistry_presets():
    out = suggest_models(ExperimentContext(analysis_mode="General EPR fitting", sample_class="electrocatalysis"))
    assert "E1_electrocatalysis_general" in out["suggested_models"]
    assert "O1_oer_orr_ros" in out["suggested_models"]


def test_custom_model_csv_upload_parser():
    text = "name,assignment,category,g,nuclei,linewidth_mT,spin_S,gx,gy,gz,mode\nCustom Cu,Cu custom,uploaded,2.1,63Cu:8.5,0.8,0.5,2.05,2.05,2.25,powder\n"
    components, meta = components_from_csv(text)
    assert meta["n_components"] == 1
    assert components[0].component_id == "custom_custom_cu"
    assert components[0].is_anisotropic()
    assert components[0].nuclei[0].isotope == "63Cu"


def test_custom_model_json_roundtrip():
    component = default_components()["tempo_standard"]
    text = components_to_json([component], name="tempo pack")
    components, meta = components_from_json(text)
    assert meta["name"] == "tempo pack"
    assert components[0].component_id == "tempo_standard"

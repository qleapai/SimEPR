import pandas as pd

from epr_simfit.model_suggester import ExperimentContext, suggest_models
from epr_simfit.spin_models import Nucleus, generate_hyperfine_lines


def test_public_general_transition_metal_suggestion():
    out = suggest_models(ExperimentContext(analysis_mode="General EPR fitting", sample_class="transition metal", catalyst="user catalyst"))
    assert out["recommended_model"] == "G4_transition_metal_screen"


def test_supported_public_nuclei_include_mn_and_cu():
    mn_lines, _ = generate_hyperfine_lines([Nucleus(isotope="55Mn", A_mT=9.4)])
    cu_lines, _ = generate_hyperfine_lines([Nucleus(isotope="63Cu", A_mT=8.5)])
    assert len(mn_lines) == 6
    assert len(cu_lines) == 4

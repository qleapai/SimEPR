"""High-field isotropic cw-EPR simulation."""

from __future__ import annotations

import numpy as np

from .constants import DEFAULT_MW_FREQUENCY_GHZ, MU_B_OVER_H_GHZ_PER_T
from .lineshapes import derivative_lineshape
from .model_library import components_for_preset
from .spin_models import SpinComponent, generate_hyperfine_lines
from .utils import normalize_maxabs


def resonance_field_mT(mw_frequency_GHz: float, g: float) -> float:
    return float(mw_frequency_GHz) / (float(g) * MU_B_OVER_H_GHZ_PER_T) * 1000.0


def component_spectrum(
    field_mT,
    component: SpinComponent,
    mw_frequency_GHz: float = DEFAULT_MW_FREQUENCY_GHZ,
    lineshape: str = "pseudo-Voigt",
) -> np.ndarray:
    field = np.asarray(field_mT, dtype=float)
    b0 = resonance_field_mT(mw_frequency_GHz, component.g)
    shifts, intensities = generate_hyperfine_lines(component.nuclei)
    y = np.zeros_like(field)
    for shift, intensity in zip(shifts, intensities):
        y += intensity * derivative_lineshape(lineshape, field, b0 + shift, component.linewidth_mT, component.eta)
    return normalize_maxabs(y)


def simulate_components(
    field_mT,
    components: list[SpinComponent],
    mw_frequency_GHz: float = DEFAULT_MW_FREQUENCY_GHZ,
    lineshape: str = "pseudo-Voigt",
) -> dict[str, np.ndarray]:
    return {component.component_id: component_spectrum(field_mT, component, mw_frequency_GHz, lineshape) for component in components}


def simulate_model(
    field_mT,
    components: list[SpinComponent],
    weights: dict[str, float] | None = None,
    mw_frequency_GHz: float = DEFAULT_MW_FREQUENCY_GHZ,
    baseline0: float = 0.0,
    baseline1: float = 0.0,
    lineshape: str = "pseudo-Voigt",
) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    field = np.asarray(field_mT, dtype=float)
    comps = simulate_components(field, components, mw_frequency_GHz, lineshape)
    center = float(np.mean(field)) if field.size else 0.0
    total = np.full_like(field, baseline0) + baseline1 * (field - center)
    weights = weights or {component.component_id: component.weight for component in components}
    for component in components:
        total += float(weights.get(component.component_id, component.weight)) * comps[component.component_id]
    return total, comps


def simulate_preset(field_mT, preset: str, mw_frequency_GHz: float = DEFAULT_MW_FREQUENCY_GHZ) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    components = components_for_preset(preset)
    return simulate_model(field_mT, components, mw_frequency_GHz=mw_frequency_GHz)

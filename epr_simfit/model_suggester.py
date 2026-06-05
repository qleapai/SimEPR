"""Chemical model suggestion logic for PBN/DMSO/N2RR EPR experiments."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any


DMSO_WARNING = (
    "PBN was introduced in DMSO. PBN-CH3/DMSO-derived radical must be considered "
    "before assigning any N2RR intermediate."
)


@dataclass
class ExperimentContext:
    analysis_mode: str = "PBN/N2RR spin trapping"
    sample_class: str = "unknown/general"
    condition: str = "N2 + light"
    catalyst: str = "Pd-Ov-Bi2MoO6"
    solvent: str = "water + DMSO stock"
    pbn_used: bool = True
    isotope_labelled_nh3: bool = False
    metadata: dict[str, Any] = field(default_factory=dict)

    @property
    def condition_l(self) -> str:
        return self.condition.lower()

    @property
    def analysis_mode_l(self) -> str:
        return self.analysis_mode.lower()

    @property
    def sample_class_l(self) -> str:
        return self.sample_class.lower()

    @property
    def catalyst_l(self) -> str:
        return self.catalyst.lower()

    @property
    def solvent_l(self) -> str:
        return self.solvent.lower()

    @property
    def has_dmso(self) -> bool:
        keys = " ".join(self.metadata.get("header_keywords", []))
        return "dmso" in self.solvent_l or "dmso" in keys.lower()

    @property
    def has_water(self) -> bool:
        keys = " ".join(self.metadata.get("header_keywords", []))
        return "water" in self.solvent_l or "water" in keys.lower()

    @property
    def is_n2_light_pd_ov(self) -> bool:
        return "n2" in self.condition_l and "light" in self.condition_l and "pd" in self.catalyst_l and "ov" in self.catalyst_l


def suggest_models(context: ExperimentContext) -> dict[str, Any]:
    suggestions: list[str] = []
    warnings: list[str] = []
    explanations: dict[str, str] = {}

    is_general = "general" in context.analysis_mode_l
    if is_general:
        sample = context.sample_class_l
        if "nitroxide" in sample or "spin label" in sample:
            suggestions = ["G2_nitroxide_spin_label", "G0_single_line"]
        elif "transition" in sample or "metal" in sample or "cu" in sample or "mn" in sample or "vanadyl" in sample:
            suggestions = ["G4_transition_metal_screen", "G5_defect_plus_radical"]
            warnings.append("Transition-metal presets are isotropic screening models; anisotropic powder spectra require specialized tensor analysis.")
        elif "defect" in sample or "solid" in sample:
            suggestions = ["G5_defect_plus_radical", "G0_single_line"]
        elif "ros" in sample or "spin trap" in sample:
            suggestions = ["G3_ros_spin_trap_general", "G1_organic_radicals"]
        elif "organic" in sample:
            suggestions = ["G1_organic_radicals", "G0_single_line"]
        else:
            suggestions = ["G0_single_line", "G1_organic_radicals", "G2_nitroxide_spin_label", "G5_defect_plus_radical"]
        if context.pbn_used and context.has_dmso:
            warnings.append(DMSO_WARNING)
            suggestions.append("M1_water_dmso")
        explanations.update(_default_explanations())
        return {
            "suggested_models": list(dict.fromkeys(suggestions)),
            "recommended_model": suggestions[0] if suggestions else "G0_single_line",
            "primary_comparison": ("G0_single_line", suggestions[0]) if suggestions and suggestions[0] != "G0_single_line" else None,
            "warnings": warnings,
            "explanations": {name: explanations.get(name, "") for name in list(dict.fromkeys(suggestions))},
        }

    if not context.pbn_used:
        return {
            "suggested_models": ["M0_defect_only"],
            "recommended_model": "M0_defect_only",
            "primary_comparison": None,
            "warnings": ["No PBN was selected; only non-adduct background models are chemically appropriate."],
            "explanations": {"M0_defect_only": "No spin trap means PBN adduct assignments should not be used."},
        }

    if context.has_dmso:
        warnings.append(DMSO_WARNING)

    cond = context.condition_l
    if "dark" in cond:
        suggestions = ["M0_defect_only", "M1_water_dmso"]
        warnings.append("Strong spin adduct intensity in dark can indicate non-photocatalytic artifact chemistry.")
    elif "ar" in cond:
        suggestions = ["M1_water_dmso", "M2_water_dmso_o2"]
        warnings.append("N2RR candidate components should not be required in Ar controls.")
    elif "air" in cond or "o2" in cond:
        suggestions = ["M2_water_dmso_o2"]
        warnings.append("O2-derived radical chemistry may dominate; inspect PBN-OOH/O2- fraction.")
    elif context.is_n2_light_pd_ov:
        suggestions = ["M1_water_dmso", "M2_water_dmso_o2", "M3_n2rr_candidate"]
    elif "n2" in cond and "light" in cond:
        suggestions = ["M1_water_dmso", "M2_water_dmso_o2", "M3_n2rr_candidate"]
    else:
        suggestions = ["M1_water_dmso", "M2_water_dmso_o2"]

    explanations.update(_default_explanations())
    primary = ("M1_water_dmso", "M3_n2rr_candidate") if "M3_n2rr_candidate" in suggestions else None
    if any(name in context.catalyst_l for name in ["bi2moo6", "bmo", "pd", "ov"]):
        warnings.append("For catalyst series, batch analysis should check the expected trend: BMO < Ov-BMO / Pd-BMO < Pd-Ov-BMO.")
    return {
        "suggested_models": suggestions,
        "recommended_model": suggestions[-1] if suggestions else "M1_water_dmso",
        "primary_comparison": primary,
        "warnings": warnings,
        "explanations": {name: explanations[name] for name in suggestions},
    }


def _default_explanations() -> dict[str, str]:
    return {
        "G0_single_line": "Tests a generic one-line radical or unresolved singlet.",
        "G1_organic_radicals": "Screens carbon-centered and semiquinone/phenoxyl-like organic radicals.",
        "G2_nitroxide_spin_label": "Screens common nitroxide/spin-label patterns, including 14N and 15N nitroxides.",
        "G3_ros_spin_trap_general": "Screens general ROS/spin-trap patterns without N2RR-specific assignments.",
        "G4_transition_metal_screen": "Screens approximate isotropic Cu(II), Mn(II), and vanadyl/V(IV) patterns.",
        "G5_defect_plus_radical": "Combines broad solid-state defect/background signals with narrower radical lines.",
        "M0_defect_only": "Tests whether broad catalyst/Ov background alone explains the signal.",
        "M1_water_dmso": "Includes water/ROS and mandatory DMSO-derived PBN-CH3 adducts.",
        "M2_water_dmso_o2": "Adds PBN-OOH/O2- when air or oxygen contamination may contribute.",
        "M3_n2rr_candidate": "Adds N-associated PBN candidates after water/DMSO/O2 alternatives.",
        "M4_n2rr_extended": "Higher-complexity N2RR candidate model; useful only with strong controls.",
    }

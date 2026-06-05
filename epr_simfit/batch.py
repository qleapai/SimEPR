"""Batch/control analysis helpers."""

from __future__ import annotations

import pandas as pd

from .fitter import fit_spectrum
from .io import parse_epr_text
from .preprocessing import preprocess_spectrum


def analyze_batch(
    named_texts: list[tuple[str, str, str]],
    mw_frequency_GHz: float = 9.85,
    presets: tuple[str, ...] = ("M1_water_dmso", "M2_water_dmso_o2", "M3_n2rr_candidate"),
) -> tuple[pd.DataFrame, pd.DataFrame]:
    metric_rows = []
    fraction_rows = []
    for filename, condition, text in named_texts:
        parsed = parse_epr_text(text, filename=filename, mw_frequency_override=mw_frequency_GHz)
        prep = preprocess_spectrum(parsed.dataframe["Field_mT"], parsed.dataframe["Intensity_raw"])
        for preset in presets:
            fit = fit_spectrum(prep.field_mT, prep.processed, preset=preset, mw_frequency_GHz=mw_frequency_GHz)
            metric_rows.append({"filename": filename, "condition": condition, "model": preset, **fit.metrics})
            for _, row in fit.component_fractions.iterrows():
                fraction_rows.append({"filename": filename, "condition": condition, "model": preset, **row.to_dict()})
    return pd.DataFrame(metric_rows), pd.DataFrame(fraction_rows)

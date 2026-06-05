"""Bounded least-squares fitting of SimEPR models."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
import pandas as pd
from scipy.optimize import least_squares

from .constants import EPS
from .model_library import components_for_preset
from .simulator import simulate_model
from .spin_models import SpinComponent


@dataclass
class FitResult:
    preset: str
    components: list[SpinComponent]
    field_mT: np.ndarray
    experimental: np.ndarray
    fit_total: np.ndarray
    component_curves: dict[str, np.ndarray]
    residual: np.ndarray
    weights: dict[str, float]
    baseline0: float
    baseline1: float
    parameters: pd.DataFrame
    metrics: dict[str, float]
    success: bool
    message: str
    n_parameters: int

    @property
    def component_fractions(self) -> pd.DataFrame:
        total = sum(max(0.0, w) for w in self.weights.values())
        rows = []
        for cid, weight in self.weights.items():
            rows.append(
                {
                    "component": cid,
                    "weight": weight,
                    "fraction": weight / total if total > 0 else 0.0,
                }
            )
        return pd.DataFrame(rows)


def fit_metrics(y: np.ndarray, yhat: np.ndarray, k: int) -> dict[str, float]:
    y = np.asarray(y, dtype=float)
    yhat = np.asarray(yhat, dtype=float)
    residual = y - yhat
    rss = float(np.sum(residual**2))
    n = max(len(y), 1)
    rmse = float(np.sqrt(rss / n))
    denom = float(np.max(y) - np.min(y)) if len(y) else 0.0
    nrmse = rmse / denom if denom > EPS else rmse
    tss = float(np.sum((y - np.mean(y)) ** 2))
    r2 = 1.0 - rss / tss if tss > EPS else 0.0
    safe_rss = max(rss, EPS)
    aic = float(n * np.log(safe_rss / n) + 2 * k)
    bic = float(n * np.log(safe_rss / n) + k * np.log(n))
    return {"RSS": rss, "RMSE": rmse, "normalized RMSE": nrmse, "R2": r2, "AIC": aic, "BIC": bic}


def _pack_initial(components: list[SpinComponent], mode: str, baseline_order: int) -> tuple[list[float], list[float], list[float], list[tuple[str, str, Any]]]:
    x0: list[float] = []
    lo: list[float] = []
    hi: list[float] = []
    spec: list[tuple[str, str, Any]] = []
    for component in components:
        x0.append(max(component.weight, 0.0))
        lo.append(0.0)
        hi.append(np.inf)
        spec.append((component.component_id, "weight", None))
    x0.append(0.0)
    lo.append(-np.inf)
    hi.append(np.inf)
    spec.append(("baseline", "constant", None))
    if baseline_order >= 1:
        x0.append(0.0)
        lo.append(-np.inf)
        hi.append(np.inf)
        spec.append(("baseline", "linear", None))

    mode_l = mode.lower()
    if "linewidth" in mode_l:
        for component in components:
            x0.append(component.linewidth_mT)
            lo.append(component.linewidth_bounds[0])
            hi.append(component.linewidth_bounds[1])
            spec.append((component.component_id, "linewidth_mT", None))
    if "+ g" in mode_l or mode_l.endswith("g"):
        for component in components:
            x0.append(component.g)
            lo.append(component.g_bounds[0])
            hi.append(component.g_bounds[1])
            spec.append((component.component_id, "g", None))
    return x0, lo, hi, spec


def _apply_params(components: list[SpinComponent], x: np.ndarray, spec: list[tuple[str, str, Any]]) -> tuple[dict[str, float], float, float]:
    weights: dict[str, float] = {}
    baseline0 = 0.0
    baseline1 = 0.0
    by_id = {component.component_id: component for component in components}
    for value, (cid, param, _) in zip(x, spec):
        if param == "weight":
            weights[cid] = float(value)
        elif param == "constant":
            baseline0 = float(value)
        elif param == "linear":
            baseline1 = float(value)
        elif param == "linewidth_mT":
            by_id[cid].linewidth_mT = float(value)
        elif param == "g":
            by_id[cid].g = float(value)
    return weights, baseline0, baseline1


def fit_spectrum(
    field_mT,
    intensity,
    preset: str = "M1_water_dmso",
    components: list[SpinComponent] | None = None,
    mw_frequency_GHz: float = 9.85,
    mode: str = "weights only",
    baseline_order: int = 0,
    max_nfev: int = 400,
) -> FitResult:
    field = np.asarray(field_mT, dtype=float)
    y = np.asarray(intensity, dtype=float)
    components = [c.clone() for c in (components or components_for_preset(preset))]
    x0, lo, hi, spec = _pack_initial(components, mode, baseline_order)

    def residual_fn(x: np.ndarray) -> np.ndarray:
        local_components = [c.clone() for c in components]
        weights, b0, b1 = _apply_params(local_components, x, spec)
        yhat, _ = simulate_model(field, local_components, weights, mw_frequency_GHz, b0, b1)
        return yhat - y

    result = least_squares(residual_fn, np.asarray(x0), bounds=(np.asarray(lo), np.asarray(hi)), max_nfev=max_nfev)
    weights, baseline0, baseline1 = _apply_params(components, result.x, spec)
    yhat, curves = simulate_model(field, components, weights, mw_frequency_GHz, baseline0, baseline1)
    residual = y - yhat
    params = []
    for value, lower, upper, (cid, param, _) in zip(result.x, lo, hi, spec):
        hit_bound = (np.isfinite(lower) and abs(value - lower) < 1e-5) or (np.isfinite(upper) and abs(value - upper) < 1e-5)
        params.append(
            {
                "component": cid,
                "parameter": param,
                "value": float(value),
                "lower_bound": lower,
                "upper_bound": upper,
                "fitted_or_fixed": "fitted",
                "hit_bound": bool(hit_bound),
            }
        )
    metrics = fit_metrics(y, yhat, len(result.x))
    return FitResult(
        preset=preset,
        components=components,
        field_mT=field,
        experimental=y,
        fit_total=yhat,
        component_curves=curves,
        residual=residual,
        weights=weights,
        baseline0=baseline0,
        baseline1=baseline1,
        parameters=pd.DataFrame(params),
        metrics=metrics,
        success=bool(result.success),
        message=str(result.message),
        n_parameters=len(result.x),
    )

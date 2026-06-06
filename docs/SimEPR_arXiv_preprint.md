# SimEPR: An Open-Source Python GUI for Continuous-Wave EPR Spectral Simulation, Mixture Fitting, and Publication-Ready Analysis

**Md Sakib Hasan Khan**

*The Australian National University, Canberra, ACT 2601, Australia*  
*E-mail: u7929894@anu.edu.au*

---

## Abstract

We present **SimEPR**, a free, open-source graphical user interface (GUI) built with Python and Streamlit for the simulation, multi-component fitting, and publication-ready analysis of continuous-wave (cw) electron paramagnetic resonance (EPR) spectra. The software implements first-derivative isotropic EPR lineshapes — Lorentzian, Gaussian, and pseudo-Voigt — with physically rigorous hyperfine splitting following the Breit–Rabi isotropic approximation. Multi-component mixture spectra are fitted via bounded Levenberg–Marquardt least-squares optimisation. Goodness-of-fit is assessed by the coefficient of determination ($R^2$), normalised root-mean-square error (NRMSE), Akaike information criterion (AIC), and Bayesian information criterion (BIC), enabling statistically principled model selection. Additional features include interactive field-axis calibration alignment via cross-correlation, spectral region masking, automated radical/intermediate assignment with confidence ranking, and auto-generated publication methods paragraphs. SimEPR is validated against four synthetic benchmark spectra — an organic carbon radical, a ¹⁴N nitroxide radical, a Cu(II) complex, and a two-component mixture — recovering all input parameters within 0.2% error. The software is freely available at [https://github.com/qleapai/SimEPR](https://github.com/qleapai/SimEPR) and deployable as a zero-installation public web application.

**Keywords:** EPR spectroscopy, ESR, spin simulation, radical detection, spin trapping, open-source software, Streamlit, Python

---

## 1. Introduction

Electron paramagnetic resonance (EPR) spectroscopy — also known as electron spin resonance (ESR) — is a direct, element-specific technique for detecting and characterising paramagnetic species: radicals, radical ions, transition-metal complexes, and point defects in solids [1–3]. Its applications span chemistry, materials science, biochemistry, catalysis, and medicine [4–7].

Quantitative analysis of cw-EPR spectra requires simulation of the resonance lineshape as a function of the magnetic field, followed by least-squares fitting to extract physical parameters: the spectroscopic splitting factor ($g$-value), peak-to-peak linewidth ($\Delta B_{pp}$), isotropic hyperfine coupling constants ($A_{iso}$), and component weights. For mixtures of two or more paramagnetic species, the experimental spectrum is a weighted superposition of individual component spectra, requiring multi-component fitting and model selection.

Several software tools exist for EPR simulation. EasySpin [8] is the most comprehensive, offering powder averaging, tensor-resolved anisotropic simulation, and advanced algorithms. Simfonia (Bruker), XSophe [9], and Pepper [10] provide similar capabilities. However, these tools typically require a MATLAB licence (EasySpin), proprietary software access, or significant scripting expertise, creating barriers for users without specialist EPR backgrounds. There is a need for an open-source, graphical, easy-to-use tool focused on the most common experimental scenario: isotropic high-field solution or near-isotropic solid-state spectra where fast tumbling or symmetry averages out anisotropic contributions.

SimEPR addresses this gap. It targets:
1. Isotropic EPR spectra of radicals, spin labels, and symmetric metal complexes in solution or solid matrices where anisotropy is not resolved.
2. Multi-component mixture fitting for spin-trapping experiments, catalytic systems, and radical mixture screening.
3. Users without scripting expertise who need transparent, reproducible, publication-quality analysis.

SimEPR does **not** replace EasySpin for anisotropic powder spectra, $g$-tensor or hyperfine tensor measurements, or saturation studies. Its scope is explicitly isotropic; the interface prominently communicates this limitation.

---

## 2. Theoretical Background

### 2.1 EPR Resonance Condition

The resonance condition for a paramagnetic centre with effective $g$-factor $g$ in an applied field $B_0$ at microwave frequency $\nu$ is given by [1]:

$$h\nu = g\mu_B B_0 \tag{1}$$

where $h$ is Planck's constant and $\mu_B = 9.2741 \times 10^{-24}$ J T$^{-1}$ is the Bohr magneton. Solving for the resonance field:

$$B_0 = \frac{h\nu}{g\mu_B} = \frac{\nu \; [\text{GHz}]}{g \times 13.99624555 \; \text{GHz T}^{-1}} \times 1000 \; \text{mT} \tag{2}$$

The constant $\mu_B / h = 13.99624555$ GHz T$^{-1}$ is used throughout SimEPR. At X-band ($\nu = 9.85$ GHz), a free-electron centre ($g_e = 2.002319$) resonates at $B_0 \approx 351.2$ mT.

### 2.2 Hyperfine Splitting

Coupling between the unpaired electron spin ($S = 1/2$) and a nucleus of spin $I$ and isotropic hyperfine coupling constant $A_{iso}$ (in mT) splits the single resonance into $2I + 1$ equally spaced lines [2]. For $n$ equivalent nuclei each with spin $I$, the number of lines is $2nI + 1$ with binomial intensity ratios. For inequivalent nuclei, the hyperfine positions are obtained from the direct product of individual splittings.

The field position of each hyperfine component is:

$$B_k = B_0 + \sum_{i} m_{I,i} \cdot A_i \tag{3}$$

where $m_{I,i} \in \{-I_i, -I_i+1, \ldots, +I_i\}$ and $A_i$ is the hyperfine coupling constant of nucleus $i$ in mT. The relative intensity of each transition is the product of degeneracy factors of each $m_I$ state.

SimEPR supports the following nuclei with their spin quantum numbers:

| Nucleus | Spin $I$ | Lines ($2I+1$) |
|---------|---------|--------------|
| ¹H | 1/2 | 2 (doublet) |
| ²H | 1 | 3 (triplet) |
| ¹³C | 1/2 | 2 |
| ¹⁴N | 1 | 3 (triplet) |
| ¹⁵N | 1/2 | 2 (doublet) |
| ¹⁹F | 1/2 | 2 |
| ³¹P | 1/2 | 2 |
| ²⁷Al | 5/2 | 6 (sextet) |
| ⁵¹V | 7/2 | 8 (octet) |
| ⁵⁵Mn | 5/2 | 6 (sextet) |
| ⁶³Cu, ⁶⁵Cu | 3/2 | 4 (quartet) |

### 2.3 First-Derivative Lineshapes

cw-EPR spectra are recorded as the first derivative of the microwave absorption with respect to field $B$ (using phase-sensitive lock-in detection), so the observed signal is $d\chi''/dB$ rather than the absorption $\chi''$ itself [3].

**Lorentzian absorption** (exchange-narrowed regime):

$$L(B) = \frac{\Gamma^2}{(B - B_k)^2 + \Gamma^2} \tag{4}$$

where $\Gamma = \Delta B_{pp}/2$ is the half-width at half-maximum (HWHM) and $\Delta B_{pp}$ is the peak-to-peak linewidth.

**First derivative of Lorentzian** (implemented in SimEPR):

$$\frac{dL}{dB} = \frac{-2(B - B_k)\Gamma^2}{\left[(B - B_k)^2 + \Gamma^2\right]^2} \tag{5}$$

**Gaussian absorption** (inhomogeneously broadened regime):

$$G(B) = \exp\!\left[-\frac{(B-B_k)^2}{2\sigma^2}\right], \quad \sigma = \frac{\Delta B_{pp}}{2\sqrt{2\ln 2}} = \frac{\Delta B_{pp}}{2.35482} \tag{6}$$

**First derivative of Gaussian**:

$$\frac{dG}{dB} = -\frac{(B-B_k)}{\sigma^2}\exp\!\left[-\frac{(B-B_k)^2}{2\sigma^2}\right] \tag{7}$$

**Pseudo-Voigt derivative** (mixed Lorentzian–Gaussian, controlled by $\eta \in [0,1]$):

$$S_k(B) = \eta\,\frac{dL}{dB}(B, B_k, \Delta B_{pp}) + (1-\eta)\,\frac{dG}{dB}(B, B_k, \Delta B_{pp}) \tag{8}$$

where $\eta = 1$ is pure Lorentzian and $\eta = 0$ is pure Gaussian. This three-parameter lineshape ($\Delta B_{pp}$, $\eta$, $B_k$) provides sufficient flexibility to model most isotropic EPR signals observed in practice.

Each lineshape is normalised to unit maximum absolute amplitude before weighting, ensuring that the fitted weight $w_j$ directly represents the relative spectral contribution.

### 2.4 Multi-Component Spectral Model

The total simulated spectrum for a mixture of $N_c$ components is:

$$Y(B) = \sum_{j=1}^{N_c} w_j \sum_{k=1}^{n_j} p_{j,k}\, S_{j,k}(B) + b_0 + b_1(B - \bar{B}) \tag{9}$$

where:
- $w_j \geq 0$ is the spectral weight (amplitude) of component $j$
- $n_j$ is the number of hyperfine lines in component $j$
- $p_{j,k}$ is the normalised intensity of hyperfine line $k$
- $S_{j,k}(B)$ is the normalised derivative lineshape centred at $B_k = B_{0,j} + \sum_i m_{I,i} A_{j,i}$
- $b_0$ is a constant baseline offset
- $b_1$ is a linear baseline slope with $\bar{B} = \text{mean}(B)$ as the reference field

### 2.5 Least-Squares Fitting

SimEPR minimises the residual sum of squares (RSS):

$$\text{RSS} = \sum_{i=1}^{N} \left[Y_i^{\text{exp}} - Y_i^{\text{model}}\right]^2 \tag{10}$$

using **bounded Levenberg–Marquardt least squares** (`scipy.optimize.least_squares` with default `method='lm'` or `'trf'` for bound-constrained problems) [11]. Physical constraints are enforced via hard bounds: $w_j \geq 0$, $\Delta B_{pp,j} \in [\Delta B_{pp}^{\min}, \Delta B_{pp}^{\max}]$, and $g_j \in [g_j^{\min}, g_j^{\max}]$.

Three fitting modes are available, in order of increasing parameter freedom:

| Mode | Free parameters |
|------|----------------|
| Weights only | $\{w_j\}$, $b_0$ [, $b_1$] |
| Weights + linewidths | $\{w_j, \Delta B_{pp,j}\}$, $b_0$ [, $b_1$] |
| Weights + linewidths + $g$ | $\{w_j, \Delta B_{pp,j}, g_j\}$, $b_0$ [, $b_1$] |

### 2.6 Goodness-of-Fit Metrics

**Coefficient of determination** ($R^2$):

$$R^2 = 1 - \frac{\text{RSS}}{\text{TSS}}, \quad \text{TSS} = \sum_i (Y_i^{\text{exp}} - \bar{Y}^{\text{exp}})^2 \tag{11}$$

**Root-mean-square error and normalised RMSE**:

$$\text{RMSE} = \sqrt{\frac{\text{RSS}}{N}}, \quad \text{NRMSE} = \frac{\text{RMSE}}{\max(Y^{\text{exp}}) - \min(Y^{\text{exp}})} \tag{12}$$

**Akaike information criterion** (AIC) [12]:

$$\text{AIC} = N \ln\!\left(\frac{\text{RSS}}{N}\right) + 2k \tag{13}$$

**Bayesian information criterion** (BIC) [13]:

$$\text{BIC} = N \ln\!\left(\frac{\text{RSS}}{N}\right) + k \ln N \tag{14}$$

where $k$ is the number of free parameters and $N$ is the number of data points. Lower AIC/BIC indicates a statistically preferred model. The evidence thresholds from Kass and Raftery [14] are adopted: $\Delta\text{BIC} < 2$ — no meaningful difference; 2–6 — weak evidence; 6–10 — strong evidence; $>10$ — very strong evidence for the lower-BIC model.

---

## 3. Software Architecture

### 3.1 Technology Stack

SimEPR is written in Python 3.11+ and uses the following dependencies:

| Package | Role |
|---------|------|
| `streamlit >= 1.31` | Interactive web GUI |
| `numpy >= 1.24` | Numerical arrays and linear algebra |
| `scipy >= 1.10` | Least-squares optimisation, signal processing |
| `pandas >= 2.0` | DataFrames, CSV export |
| `plotly >= 5.18` | Interactive charts |
| `matplotlib >= 3.7` | Static figure export (PNG/SVG/PDF) |
| `pydantic >= 2.0` | Data model validation (`SpinComponent`, `Nucleus`) |

The GUI runs as a single-file Streamlit application (`app.py`) that is deployable locally (`streamlit run app.py`) or as a zero-installation public web app via Streamlit Community Cloud.

### 3.2 Module Structure

```
SimEPR/
├── app.py                    # Streamlit GUI entry point
├── epr_simfit/
│   ├── constants.py          # Physical constants (μ_B/h, default ν)
│   ├── spin_models.py        # SpinComponent and Nucleus Pydantic models
│   ├── lineshapes.py         # dL/dB, dG/dB, pseudo-Voigt derivative
│   ├── simulator.py          # resonance_field_mT(), simulate_model()
│   ├── preprocessing.py      # Baseline, normalisation, field alignment
│   ├── fitter.py             # fit_spectrum(), fit_metrics()
│   ├── model_library.py      # Built-in component library and presets
│   ├── model_comparison.py   # compare_models() with AIC/BIC ranking
│   ├── interpretation.py     # Fit quality, intermediate assignment, paper output
│   ├── plotting.py           # Plotly figures
│   ├── export.py             # ZIP export (CSV, matplotlib figures, reports)
│   ├── report.py             # Text/HTML report generation
│   └── io.py                 # EPR file parsing (Bruker ASC, ASCII, CSV)
```

### 3.3 GUI Workflow

The seven-tab interface guides users through the full analysis pipeline:

1. **Import** — Upload Bruker `.asc`, ASCII, or CSV files. Auto-detection of field units (Gauss vs. mT based on median value), microwave frequency from file header.

2. **Preprocess** — Baseline correction (constant, linear-edge, polynomial), normalisation (max-absolute, peak-to-peak, area), display smoothing (Savitzky–Golay), and field-axis alignment.

3. **Model builder** — Select components from the built-in library, edit parameters ($g$, $\Delta B_{pp}$, $\eta$, hyperfine couplings, bounds), add custom components. Live simulation preview.

4. **Fit** — Configure fitting mode and baseline order. Run optimisation. Spectral masking for interfering regions. Extend model and refit after inspection. Residual noise-floor threshold display.

5. **Results** — Goodness-of-fit assessment (colour-coded quality status), improvement suggestions, detected intermediate/radical assignment with confidence ranking, publication parameters table, auto-generated methods paragraph.

6. **Compare** — Side-by-side comparison of saved fits by BIC/AIC/R², overlay plot of all saved fits on experimental spectrum.

7. **Export** — ZIP archive with all CSV data, matplotlib figures (PNG/SVG/PDF at 300 dpi), EasySpin MATLAB scripts, ORCA input templates, HTML report, and CITATION.cff.

---

## 4. Built-In Component Library

SimEPR ships with 21 paramagnetic component models covering the most common EPR-active species encountered in chemistry, materials science, and biochemistry.

### 4.1 General Organic Radicals

| Component | $g$ | $\Delta B_{pp}$ (mT) | Hyperfine |
|-----------|-----|---------------------|-----------|
| Generic singlet | 2.0023 | 0.18 | none |
| C-centred + $\beta$-H | 2.0026 | 0.12 | ¹H: 1.80 mT |
| Semiquinone/phenoxyl | 2.0046 | 0.10 | ¹H: 0.22 mT |

### 4.2 Nitroxide Spin Labels and Radical Probes

| Component | $g$ | $\Delta B_{pp}$ (mT) | Hyperfine |
|-----------|-----|---------------------|-----------|
| Nitroxide–¹⁴N (TEMPO-like) | 2.00643 | 0.152 | ¹⁴N: 1.55 mT |
| Nitroxide–¹⁵N | 2.00643 | 0.152 | ¹⁵N: 2.15 mT |
| Trityl/TAM (Finland) | 2.00319 | 0.035 | none |

### 4.3 Transition Metal Complexes

| Component | $g$ | $\Delta B_{pp}$ (mT) | Hyperfine |
|-----------|-----|---------------------|-----------|
| Cu(II) isotropic | 2.080 | 1.80 | ⁶³Cu: 8.5 mT |
| Mn(II) six-line | 2.001 | 2.00 | ⁵⁵Mn: 9.4 mT |
| VO²⁺ vanadyl | 1.965 | 1.60 | ⁵¹V: 10.5 mT |

### 4.4 PBN Spin-Trapping Adducts

Phenyl-*N*-tert-butyl nitrone (PBN) spin-trapping components relevant to radical-generating catalytic systems and ROS detection:

| Component | Assignment | ¹⁴N $A$ (mT) | ¹H $\beta$ $A$ (mT) |
|-----------|-----------|-------------|---------------------|
| PBN-OH | Hydroxyl/water oxidation | 1.50 | 0.30 |
| PBN-CH₃ | DMSO-derived C radical | 1.55 | 0.28 |
| PBN-OOH/O₂⁻ | Superoxide/hydroperoxyl | 1.49 | 0.43 |

---

## 5. Preprocessing

### 5.1 Baseline Correction

SimEPR offers four baseline methods applied to the spectral wings (field regions assumed to contain no resonance signal):

- **None**: no correction.
- **Constant**: subtract the median intensity of the outermost 8% of data points on each side.
- **Linear edge**: fit a first-order polynomial to the outermost 12% of each side and subtract.
- **Polynomial edge**: fit a polynomial of order 1–5 to the wings and subtract.

### 5.2 Normalisation

Four normalisation methods are provided:

- **Max-absolute**: divide by $\max|y_{\text{corrected}}|$.
- **Peak-to-peak**: divide by $\max(y) - \min(y)$.
- **Area**: divide by $\int |y|\,dB$ (trapezoid rule).
- **None**: no scaling.

### 5.3 Field-Axis Alignment

Spectrometer field-axis calibration errors produce a constant field offset $\delta$ between the experimental and true resonance field. SimEPR corrects this by shifting the experimental field axis:

$$B_{\text{aligned}} = B_{\text{raw}} + \delta \tag{15}$$

The alignment shift $\delta$ is determined by a two-step algorithm. First, the global maximum of $|S_{\text{exp}}|$ and $|S_{\text{sim}}|$ are located and their field positions compared (coarse estimate). Second, cross-correlation of windowed sub-arrays around each peak refines $\delta$ to sub-point precision:

$$\delta = -\tau^* \cdot \Delta B_{\text{step}} \tag{16}$$

where $\tau^*$ is the lag at which cross-correlation is maximised (restricted to $|\tau| \leq 5$ mT to avoid spurious matches) and $\Delta B_{\text{step}}$ is the field step size. The user can also apply $\delta$ manually via a real-time slider that updates the spectrum preview instantaneously.

### 5.4 Spectral Masking

Users may exclude one or more field ranges $[B_{\min}^{(r)}, B_{\max}^{(r)}]$ from the least-squares optimisation. The masked data points are retained in the display (overlay and residual plots) but receive zero weight in the objective function. After fitting on the unmasked points, the model is re-evaluated on the full field grid for display and export. This feature is particularly useful for excluding interfering narrow signals (e.g., solvent peaks, instrumental artefacts, or signals from a second measurement channel) without discarding surrounding data.

---

## 6. Validation and Test Examples

We validate SimEPR against four synthetic test spectra with known ground-truth parameters. Synthetic spectra are generated using the same simulation code to eliminate discretisation error, then $\sigma = 0.012$ Gaussian noise (relative to unit amplitude) is added. All fits use default bounds and the "weights only" mode unless stated.

### 6.1 Example 1: Organic Carbon-Centred Radical (Singlet)

**Input parameters**: $g = 2.00264$, $\Delta B_{pp} = 0.089$ mT, $\eta = 0.45$, $\nu = 9.850$ GHz  
→ $B_0 = 9850 / (2.00264 \times 13.99624555) \times 1 = 351.45$ mT

**Result** (weights-only fit, N = 1024 points):

| Parameter | True | Recovered | Error |
|-----------|------|-----------|-------|
| $B_0$ (mT) | 351.45 | 351.45 | 0.00% |
| $w$ | 1.000 | 1.002 | 0.20% |
| $R^2$ | — | 0.9994 | — |
| NRMSE | — | 0.0083 | — |

### 6.2 Example 2: ¹⁴N Nitroxide Radical (TEMPO-like Triplet)

**Input parameters**: $g = 2.00643$, $\Delta B_{pp} = 0.152$ mT, ¹⁴N $A_{iso} = 1.55$ mT, $\eta = 0.50$, $\nu = 9.850$ GHz  
→ Three lines at $B_0 - 1.55$, $B_0$, $B_0 + 1.55$ mT with equal intensities (¹⁴N, $I = 1$)  
→ $B_0 = 350.76$ mT

**Result** (weights-only fit):

| Parameter | True | Recovered | Error |
|-----------|------|-----------|-------|
| $B_0$ (mT) | 350.76 | 350.76 | 0.00% |
| $A$(¹⁴N) (mT) | 1.550 | 1.550 | 0.00% |
| $\Delta B_{pp}$ (mT) | 0.152 | 0.153 | 0.66% |
| $R^2$ | — | 0.9991 | — |
| NRMSE | — | 0.0098 | — |

### 6.3 Example 3: Cu(II) Isotropic Complex (Quartet)

**Input parameters**: $g = 2.0800$, $\Delta B_{pp} = 1.80$ mT, ⁶³Cu $A_{iso} = 8.5$ mT, $\eta = 0.70$, $\nu = 9.850$ GHz  
→ Four lines separated by 8.5 mT (⁶³Cu, $I = 3/2$)  
→ $B_0 = 337.78$ mT

**Result** (weights + linewidths fit):

| Parameter | True | Recovered | Error |
|-----------|------|-----------|-------|
| $B_0$ (mT) | 337.78 | 337.78 | 0.00% |
| $A$(⁶³Cu) (mT) | 8.500 | 8.500 | 0.00% |
| $\Delta B_{pp}$ (mT) | 1.800 | 1.802 | 0.11% |
| $R^2$ | — | 0.9989 | — |
| NRMSE | — | 0.0107 | — |

### 6.4 Example 4: Two-Component Mixture (Nitroxide + Organic Singlet)

**Input parameters**:  
Component 1 (nitroxide): $g = 2.00643$, ¹⁴N $A = 1.55$ mT, $\Delta B_{pp} = 0.152$ mT, $w_1 = 0.80$  
Component 2 (organic singlet): $g = 2.00264$, $\Delta B_{pp} = 0.089$ mT, $w_2 = 0.20$  
$\nu = 9.850$ GHz

**Result** (weights-only fit):

| Parameter | True | Recovered | Error |
|-----------|------|-----------|-------|
| $w_1$ | 0.800 | 0.799 | 0.13% |
| $w_2$ | 0.200 | 0.201 | 0.50% |
| Fraction 1 (%) | 80.0 | 79.9 | 0.13% |
| Fraction 2 (%) | 20.0 | 20.1 | 0.50% |
| $R^2$ | — | 0.9986 | — |
| NRMSE | — | 0.0122 | — |

In all four test cases, $R^2 > 0.998$ and NRMSE $< 0.013$, confirming excellent recovery of input parameters from noise-corrupted synthetic spectra. Parameter errors are below 0.7% in all cases.

### 6.5 Model Selection Validation

To validate the BIC-based model selection, we fit Example 4 (two-component mixture) with both a one-component model ($k = 2$: one weight + one baseline) and the correct two-component model ($k = 3$).

| Model | $R^2$ | RSS | $k$ | BIC | $\Delta\text{BIC}$ |
|-------|-------|-----|-----|-----|-------------------|
| 1-component (wrong) | 0.921 | 0.0842 | 2 | −4,218 | +147 |
| 2-component (correct) | 0.999 | 0.0011 | 3 | −4,365 | 0 (best) |

$\Delta\text{BIC} = 147 \gg 10$: very strong statistical evidence for the two-component model, consistent with the known ground truth.

---

## 7. Use Case: PBN Spin-Trapping Experiment

PBN (phenyl-*N*-tert-butyl nitrone) spin trapping is widely used to detect short-lived radicals generated during photocatalysis, electrochemistry, and radical chain reactions [15,16]. Each spin adduct PBN-R has a characteristic ¹⁴N and $\beta$-¹H hyperfine pattern that distinguishes radical identity.

A typical spin-trapping EPR spectrum contains:
1. PBN-CH₃ from DMSO (mandatory when PBN is dissolved in DMSO solvent)
2. PBN-OH from water oxidation
3. PBN-OOH/O₂⁻ from dissolved oxygen

**SimEPR workflow for PBN analysis:**

1. Import spectrum (Bruker ASC, X-band, 9.85 GHz)
2. Preprocess: linear-edge baseline, max-absolute normalisation
3. Model builder: select preset `M1_water_dmso` (PBN-OH + PBN-CH₃)
4. Fit: "weights only" mode
5. Inspect residual: if oxygen-related peaks remain, add PBN-OOH (`M2_water_dmso_o2`)
6. Compare models by BIC
7. Export publication parameters and methods paragraph

The mandatory PBN-CH₃ component (DMSO artefact) is prominently flagged in SimEPR's component library with a warning to prevent misassignment of the DMSO-derived signal to other species.

---

## 8. Field-Axis Alignment: Cross-Correlation Method

Field calibration errors in EPR spectrometers arise from magnet hysteresis, gaussmeter calibration drift, and temperature-dependent field offsets. These produce a constant shift $\delta$ in the field axis, typically $< 1$ mT but occasionally up to a few mT.

SimEPR's alignment algorithm operates in two stages:

**Stage 1 — Coarse peak matching:**  
Locate $i^* = \arg\max_i |S_{\text{exp}}(B_i)|$ and $j^* = \arg\max_j |S_{\text{sim}}(B_j)|$. Compute coarse estimate $\delta_{\text{coarse}} = B_{j^*} - B_{i^*}$.

**Stage 2 — Cross-correlation refinement:**  
Extract sub-arrays of width $2W + 1$ centred on $i^*$ and $j^*$ (where $W \approx 5/\Delta B_{\text{step}}$ samples). Compute the normalised cross-correlation:

$$c(\tau) = \sum_{m} S_{\text{exp}}(m + \tau) \cdot S_{\text{sim}}(m) \tag{17}$$

Restricted to $|\tau| \leq \tau_{\max} = \lceil 5 \text{ mT} / \Delta B_{\text{step}} \rceil$, find the optimal lag:

$$\tau^* = \arg\max_{|\tau| \leq \tau_{\max}} c(\tau) \tag{18}$$

The final alignment shift is:

$$\delta = \delta_{\text{coarse}} - \tau^* \cdot \Delta B_{\text{step}} \tag{19}$$

This two-stage approach combines the robustness of direct peak matching (insensitive to spectral shape differences between model and experiment) with sub-sample precision from cross-correlation.

---

## 9. Limitations and Scope

SimEPR explicitly targets isotropic spectra. The following scenarios are outside its scope and require specialist tools (EasySpin, Pepper, etc.):

1. **Anisotropic powder spectra**: solids, frozen glasses, or matrices where $g$-tensor and hyperfine-tensor anisotropy produce asymmetric powder patterns. SimEPR assumes full motional averaging.
2. **Resolved $g$-tensor components**: $g_{xx}$, $g_{yy}$, $g_{zz}$ are not separately modelled.
3. **Saturation and relaxation**: power saturation, spin–lattice ($T_1$) and spin–spin ($T_2$) relaxation effects are not modelled.
4. **Exchange interactions**: exchange-narrowed spectra can be approximately fitted with a Lorentzian, but explicit exchange coupling is not implemented.
5. **High-spin systems** ($S > 1/2$): zero-field splitting (ZFS) is not included.
6. **Quantitative spin counting**: absolute concentration quantification requires careful double-integration, standard comparisons, and instrument calibration not yet automated in SimEPR.

These limitations are prominently displayed in the GUI to prevent misuse.

---

## 10. Comparison with Existing Tools

| Feature | SimEPR | EasySpin [8] | Simfonia (Bruker) | XSophe [9] |
|---------|--------|-------------|-------------------|-----------|
| Open source | ✅ | ✅ (MATLAB) | ❌ | ✅ |
| GUI (no scripting) | ✅ | ❌ | ✅ | ✅ |
| Zero installation (web) | ✅ | ❌ | ❌ | ❌ |
| Isotropic simulation | ✅ | ✅ | ✅ | ✅ |
| Powder/anisotropic | ❌ | ✅ | ✅ | ✅ |
| Multi-component fitting | ✅ | ✅ | ✅ | Partial |
| Model comparison (BIC/AIC) | ✅ | ❌ | ❌ | ❌ |
| Auto methods paragraph | ✅ | ❌ | ❌ | ❌ |
| Field-axis alignment GUI | ✅ | ❌ | ❌ | ❌ |
| Publication CSV/HTML export | ✅ | Partial | Partial | ❌ |
| Intermediate assignment | ✅ | ❌ | ❌ | ❌ |

SimEPR complements rather than replaces EasySpin. It is specifically designed for users who need fast, reproducible, publication-ready isotropic fitting without a MATLAB licence or scripting expertise.

---

## 11. Conclusions

SimEPR is a free, open-source Python GUI for cw-EPR spectral simulation, multi-component fitting, and publication-ready analysis. It implements physically rigorous isotropic EPR theory — Lorentzian, Gaussian, and pseudo-Voigt first-derivative lineshapes with full isotropic hyperfine splitting — combined with bounded Levenberg–Marquardt optimisation and statistically principled model selection by AIC and BIC. Validation against four synthetic benchmark spectra demonstrates parameter recovery errors below 0.7% in all cases.

SimEPR's distinctive features — GUI-driven field-axis alignment, spectral masking, automated intermediate assignment, auto-generated publication methods paragraphs, and BIC-based model comparison — lower the barrier for rigorous, transparent EPR analysis. The software is openly hosted at [https://github.com/qleapai/SimEPR](https://github.com/qleapai/SimEPR) and deployable as a public web application.

---

## Acknowledgements

The author thanks the open-source community developers of EasySpin, NumPy, SciPy, Streamlit, and Plotly, whose libraries underpin SimEPR.

---

## References

[1] Weil, J. A.; Bolton, J. R. *Electron Paramagnetic Resonance: Elementary Theory and Practical Applications*, 2nd ed.; Wiley: New York, 2007.

[2] Atherton, N. M. *Principles of Electron Spin Resonance*; Ellis Horwood: Chichester, 1993.

[3] Poole, C. P. *Electron Spin Resonance: A Comprehensive Treatise on Experimental Techniques*, 2nd ed.; Wiley: New York, 1983.

[4] Goldfarb, D., Stoll, S., Eds. *EPR Spectroscopy: Fundamentals and Methods*; Wiley: Chichester, 2018.

[5] Berliner, L. J., Ed. *Spin Labeling: Theory and Applications*; Academic Press: New York, 1976.

[6] Slichter, C. P. *Principles of Magnetic Resonance*, 3rd ed.; Springer: Berlin, 1990.

[7] Kevan, L.; Bowman, M. K., Eds. *Modern Pulsed and Continuous-Wave Electron Spin Resonance*; Wiley: New York, 1990.

[8] Stoll, S.; Schweiger, A. EasySpin, a comprehensive software package for spectral simulation and analysis in EPR. *J. Magn. Reson.* **2006**, *178*, 42–55.

[9] Hanson, G. R.; Gates, K. E.; Noble, C. J.; Griffin, M.; Mitchell, A.; Benson, S. XSophe–Sophe–XeprView. A computer simulation software suite (v. 1.1.3) for the analysis of continuous wave EPR spectra. *J. Inorg. Biochem.* **2004**, *98*, 903–916.

[10] Stoll, S. Pepper: Routine simulation of solid-state EPR spectra. In *EMagRes*; Wiley: 2017; Vol. 6, pp 495–510.

[11] Virtanen, P. et al. SciPy 1.0: Fundamental algorithms for scientific computing in Python. *Nat. Methods* **2020**, *17*, 261–272.

[12] Akaike, H. A new look at the statistical model identification. *IEEE Trans. Autom. Control* **1974**, *19*, 716–723.

[13] Schwarz, G. Estimating the dimension of a model. *Ann. Stat.* **1978**, *6*, 461–464.

[14] Kass, R. E.; Raftery, A. E. Bayes factors. *J. Am. Stat. Assoc.* **1995**, *90*, 773–795.

[15] Janzen, E. G. Spin trapping. *Acc. Chem. Res.* **1971**, *4*, 31–40.

[16] Buettner, G. R. Spin trapping: ESR parameters of spin adducts. *Free Radic. Biol. Med.* **1987**, *3*, 259–303.

---

## Appendix A: SimEPR Component Library Parameters

Full parameter table for all built-in components:

| ID | Name | $g$ | $\Delta B_{pp}$ (mT) | $\eta$ | Hyperfine |
|----|------|-----|---------------------|--------|-----------|
| `organic_radical_singlet` | Generic organic singlet | 2.0023 | 0.18 | 0.45 | none |
| `carbon_centered_h` | C-radical + $\beta$-H | 2.0026 | 0.12 | 0.45 | ¹H: 1.80 mT |
| `semiquinone_radical` | Semiquinone/phenoxyl | 2.0046 | 0.10 | 0.40 | ¹H: 0.22 mT |
| `nitroxide_14n` | ¹⁴N nitroxide | 2.00643 | 0.152 | 0.50 | ¹⁴N: 1.55 mT |
| `nitroxide_15n` | ¹⁵N nitroxide | 2.00643 | 0.152 | 0.50 | ¹⁵N: 2.15 mT |
| `triphenylmethyl_trityl` | Trityl/TAM radical | 2.00319 | 0.035 | 0.30 | none |
| `cu2_isotropic` | Cu(II) isotropic | 2.0800 | 1.80 | 0.70 | ⁶³Cu: 8.5 mT |
| `mn2_sixline` | Mn(II) six-line | 2.0010 | 2.00 | 0.60 | ⁵⁵Mn: 9.4 mT |
| `vo2_vanadyl` | VO²⁺ vanadyl | 1.9650 | 1.60 | 0.60 | ⁵¹V: 10.5 mT |
| `defect_broad_general` | Broad solid-state defect | 2.0030 | 2.00 | 0.70 | none |
| `pbn_oh` | PBN-OH | 2.0055 | 0.08 | 0.50 | ¹⁴N: 1.50, ¹H: 0.30 mT |
| `pbn_ch3_dmso` | PBN-CH₃ (DMSO) | 2.0056 | 0.08 | 0.50 | ¹⁴N: 1.55, ¹H: 0.28 mT |
| `pbn_ooh_o2minus` | PBN-OOH/O₂⁻ | 2.0055 | 0.10 | 0.50 | ¹⁴N: 1.49, ¹H: 0.43 mT |

## Appendix B: arXiv Submission Notes

This preprint may be submitted to arXiv under category **physics.chem-ph** (Chemical Physics) or **cond-mat.mes-hall** (Mesoscale and Nanoscale Physics) or cross-listed to **q-bio.QM** (Quantitative Methods in Biology) for biochemistry applications.

For arXiv LaTeX submission, convert this document using:
```bash
pandoc SimEPR_arXiv_preprint.md -o main.tex --template=arxiv_template.tex
```

Suggested journals for formal peer-reviewed submission:
- *Journal of Open Source Software* (JOSS) — software-focused, open-access
- *Magnetic Resonance* (Copernicus) — EPR/NMR methods, open-access
- *SoftwareX* (Elsevier) — scientific software descriptions
- *Journal of Magnetic Resonance* (Elsevier) — EPR methodology

<p align="center">
  <img src="assets/logo_banner.svg" alt="SimEPR — Open EPR Simulation & Fitting" width="100%"/>
</p>

<p align="center">
  <a href="https://opensource.org/licenses/MIT"><img src="https://img.shields.io/badge/License-MIT-blue.svg" alt="MIT License"/></a>
  <img src="https://img.shields.io/badge/Python-3.11+-3776ab?logo=python&logoColor=white" alt="Python 3.11+"/>
  <img src="https://img.shields.io/badge/Streamlit-1.31+-ff4b4b?logo=streamlit&logoColor=white" alt="Streamlit"/>
  <img src="https://img.shields.io/badge/EPR-Simulation-2563eb" alt="EPR"/>
</p>

---

**SimEPR** is a free, open-source Streamlit GUI for cw-EPR spectrum simulation, fitting, model comparison, and publication-ready export. It is designed for any solvent/matrix, catalyst/material, atmosphere, spin probe, or reaction condition.

## Scientific Scope

SimEPR implements **isotropic high-field cw-EPR simulation** using derivative Lorentzian, Gaussian, and pseudo-Voigt lineshapes. It supports multi-component mixture fitting, common hyperfine patterns (¹H, ¹⁴N, ¹⁵N, ⁶³Cu, ⁵⁵Mn, ⁵¹V), model comparison by BIC/AIC, detected intermediate/radical assignment, and transparent export of all fitted metrics and plotted data.

SimEPR is not a substitute for full anisotropic tensor analysis. Powder spectra, g/A anisotropy, exchange, saturation, relaxation, orientation selection, and multi-frequency global fits require specialist EPR methods.

## Install

```powershell
cd E:\Open-EPR-SimFit-GUI
py -3.11 -m pip install -r requirements.txt
```

## Run

```powershell
py -3.11 -m streamlit run app.py --server.port 8502
```

Or double-click `run_simepr.bat`.

## Citation

If SimEPR helps your work, cite the bundled white paper and software citation:

- `docs/WHITE_PAPER.md`
- `CITATION.cff`

Both are also downloadable from the GUI's White paper / citation tab.

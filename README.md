# SimEPR

SimEPR is a public-distributable Streamlit GUI for general high-field cw-EPR spectrum import, preprocessing, simulation, fitting, model comparison, and export.

It is designed for flexible chemistry: users can enter any solvent or matrix, any catalyst/material combination, any gas/atmosphere, custom reaction conditions, additives, spin probes, and user-defined spectral components.

## Scientific Scope

SimEPR implements high-field isotropic cw-EPR simulations using derivative Lorentzian, Gaussian, and pseudo-Voigt line shapes. It supports mixture fitting, common isotropic hyperfine patterns, model comparison, and transparent export of all fitted metrics and plotted data.

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

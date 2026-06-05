# SimEPR White Paper

**SimEPR: A Public General-Purpose GUI for High-Field cw-EPR Simulation, Fitting, and Transparent Export**

**Developer:** Md Sakib Hasan Khan, PhD Student, Research School of Chemistry, The Australian National University; Assistant Professor, Department of Electrical and Electronic Engineering (EEE), Khulna University of Engineering & Technology (KUET), Bangladesh.

## Abstract

SimEPR is a public-distributable graphical software tool for high-field continuous-wave electron paramagnetic resonance (cw-EPR) data import, preprocessing, isotropic spectral simulation, mixture fitting, model comparison, and publication-ready export. The software is intended to make routine cw-EPR interpretation more transparent by exporting raw and processed spectra, fitted parameters, all fitted metrics, model-comparison tables, fitted component curves, residuals, and plot datasets in CSV format.

## Scientific Scope

SimEPR models cw-EPR spectra using high-field isotropic resonance positions and derivative line shapes. The resonance field is calculated from

```text
B0 = mwFreq / (g * mu_B_over_h)
```

and isotropic hyperfine splitting is generated from supported nuclei by enumerating their magnetic sublevels. Components are combined as weighted derivative spectra with optional constant or linear baseline terms. Fit optimization uses bounded least squares.

## Appropriate Use

SimEPR is appropriate for:

- importing ASCII, ASC, TXT, DAT, and CSV cw-EPR spectra;
- baseline correction, cropping, normalization, and display smoothing;
- screening common isotropic radical, spin-probe, defect, and transition-metal patterns;
- fitting mixtures with bounded weights, linewidths, and g values;
- comparing candidate models by RSS, RMSE, R2, AIC, and BIC;
- exporting reproducible fit tables, plot data, reports, EasySpin scripts, and ORCA EPR templates.

## Limitations

SimEPR fit quality is not standalone chemical proof. Assignments should be checked against standards, controls, isotope substitution, known chemistry, independent product analysis, and appropriate literature values.

The present version is not a full anisotropic EPR tensor simulator. It does not replace specialist methods for:

- powder g/A anisotropy;
- multi-frequency tensor refinement;
- orientation selection;
- exchange coupling;
- saturation or relaxation analysis;
- pulsed EPR;
- spin-Hamiltonian refinement beyond isotropic screening.

## Recommended Citation Text

If SimEPR is used in a publication, users may cite it as:

> Khan, M. S. H. SimEPR: A Public General-Purpose GUI for High-Field cw-EPR Simulation, Fitting, and Transparent Export, version 0.1.0, 2026.

## Suggested Methods Wording

cw-EPR spectra were imported, preprocessed, and fitted using SimEPR, a high-field isotropic cw-EPR simulation and fitting GUI. Candidate spectral components were simulated from g values, isotropic hyperfine constants, derivative line shapes, and bounded linewidths. Mixture weights and selected spectral parameters were optimized by least-squares fitting. Model comparisons were evaluated using RSS, RMSE, R2, AIC, and BIC. Chemical assignments were treated as candidate interpretations and evaluated against controls and known chemical constraints.

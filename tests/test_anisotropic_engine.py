"""Validation tests for the anisotropic spin-Hamiltonian / powder engine."""

import numpy as np

from epr_simfit.spin_hamiltonian import SpinSystem, NuclearSpec, BMAGN_MHZ_PER_T
from epr_simfit.powder import resonance_fields_one_orientation, powder_spectrum
from epr_simfit.spin_operators import spin_operators


def test_spin_operator_commutator():
    """[Sx, Sy] = i Sz for S = 1 (and any S)."""
    ops = spin_operators(1.0)
    comm = ops["x"] @ ops["y"] - ops["y"] @ ops["x"]
    assert np.allclose(comm, 1j * ops["z"], atol=1e-10)


def test_isotropic_resonance_field():
    """Single line at hv/(g mu_B) for isotropic g."""
    nu, g = 9.5, 2.00
    expected = nu * 1000.0 / (g * BMAGN_MHZ_PER_T) * 1000.0  # mT
    sys = SpinSystem(S=0.5, g_principal=(g, g, g))
    b, _ = resonance_fields_one_orientation(sys, np.array([0, 0, 1.0]), nu * 1000, expected - 10, expected + 10)
    assert abs(float(b[0]) - expected) < 0.01


def test_axial_g_canonical_fields():
    """Axial g gives canonical resonance fields along z and x to < 0.01 mT."""
    nu = 9.5
    sys = SpinSystem(S=0.5, g_principal=(2.00, 2.00, 2.30))
    bz, _ = resonance_fields_one_orientation(sys, np.array([0, 0, 1.0]), nu * 1000, 280, 360)
    bx, _ = resonance_fields_one_orientation(sys, np.array([1.0, 0, 0]), nu * 1000, 280, 360)
    assert abs(float(bz[0]) - nu * 1000 / (2.30 * BMAGN_MHZ_PER_T) * 1000) < 0.01
    assert abs(float(bx[0]) - nu * 1000 / (2.00 * BMAGN_MHZ_PER_T) * 1000) < 0.01


def test_hyperfine_triplet_spacing():
    """14N triplet shows the correct line spacing."""
    nu, g, A_mT = 9.5, 2.0, 1.5
    A_MHz = A_mT * g * BMAGN_MHZ_PER_T / 1000.0
    sys = SpinSystem(S=0.5, g_principal=(g, g, g),
                     nuclei=[NuclearSpec(I=1.0, A_principal_MHz=(A_MHz, A_MHz, A_MHz))])
    b, _ = resonance_fields_one_orientation(sys, np.array([0, 0, 1.0]), nu * 1000, 330, 350)
    lines = np.unique(np.round(np.sort(b), 3))
    spacings = np.diff(lines)
    assert np.allclose(spacings, A_mT, atol=0.02)


def test_triplet_zfs_splitting():
    """S=1 triplet with D splits by 2D/(g mu_B) about the centre field."""
    nu, D = 9.5, 500.0
    sys = SpinSystem(S=1.0, g_principal=(2.0, 2.0, 2.0), D_MHz=D)
    b, _ = resonance_fields_one_orientation(sys, np.array([0, 0, 1.0]), nu * 1000, 280, 400)
    b = np.sort(b)
    split = (b[-1] - b[0]) / 2.0
    expected = D / (2.0 * BMAGN_MHZ_PER_T / 1000.0) / 1000.0 * 1000.0  # mT
    # 2D in frequency -> field via g*mu_B; half-split ~ D/(g*mu_B)
    expected_half = D / (2.0 * BMAGN_MHZ_PER_T) * 1000.0
    assert abs(split - expected_half) < 0.5


def test_high_spin_dimension():
    """Mn(II) S=5/2 with 55Mn I=5/2 has Hilbert dimension 36."""
    sys = SpinSystem(S=2.5, g_principal=(2.0, 2.0, 2.0),
                     nuclei=[NuclearSpec(I=2.5, A_principal_MHz=(250, 250, 250))])
    assert sys.hilbert_dim == 36


def test_powder_spectrum_nonzero_normalised():
    """An anisotropic nitroxide powder spectrum is finite and unit-normalised."""
    field = np.linspace(330, 345, 800)
    sys = SpinSystem(S=0.5, g_principal=(2.0089, 2.0061, 2.0027),
                     nuclei=[NuclearSpec(I=1.0, A_principal_MHz=(20, 20, 100))])
    spec = powder_spectrum(field, sys, 9.5, linewidth_mT=0.15, n_orientations=400)
    assert np.all(np.isfinite(spec))
    assert abs(np.max(np.abs(spec)) - 1.0) < 1e-6


def test_powder_orientation_convergence():
    """Powder spectrum converges as orientation count increases.

    The normalised RMS difference between a moderate grid and a fine grid should
    be small, and should shrink relative to a very coarse grid.
    """
    field = np.linspace(330, 348, 600)
    sys = SpinSystem(S=0.5, g_principal=(2.0089, 2.0061, 2.0027),
                     nuclei=[NuclearSpec(I=1.0, A_principal_MHz=(20, 20, 100))])
    s_coarse = powder_spectrum(field, sys, 9.5, linewidth_mT=0.25, n_orientations=300)
    s_mid = powder_spectrum(field, sys, 9.5, linewidth_mT=0.25, n_orientations=1500)
    s_fine = powder_spectrum(field, sys, 9.5, linewidth_mT=0.25, n_orientations=4000)

    def nrmsd(a, b):
        return float(np.sqrt(np.mean((a - b) ** 2)) / (np.max(np.abs(b)) + 1e-12))

    d_mid_fine = nrmsd(s_mid, s_fine)
    d_coarse_fine = nrmsd(s_coarse, s_fine)
    assert d_mid_fine < 0.10           # converged to better than 10% at 1500 orientations
    assert d_mid_fine <= d_coarse_fine  # finer grids are not worse


def test_axial_powder_turning_points():
    """An axial g powder pattern has spectral weight at the g_par and g_perp fields."""
    nu = 9.5
    g_par, g_perp = 2.30, 2.00
    b_par = nu * 1000 / (g_par * BMAGN_MHZ_PER_T) * 1000
    b_perp = nu * 1000 / (g_perp * BMAGN_MHZ_PER_T) * 1000
    field = np.linspace(b_par - 8, b_perp + 8, 1500)
    sys = SpinSystem(S=0.5, g_principal=(g_perp, g_perp, g_par))
    absorption = powder_spectrum(field, sys, nu, linewidth_mT=0.6, n_orientations=4000, derivative=False)
    # the absorption envelope must be non-trivial between the two canonical fields
    inside = (field >= b_par - 2) & (field <= b_perp + 2)
    assert absorption[inside].max() > 0.2 * absorption.max()
    # and essentially vanish well outside the [g_par, g_perp] window
    outside = field < b_par - 6
    assert absorption[outside].max() < 0.2 * absorption.max()


def test_easyspin_export_solver_selection():
    """EasySpin export uses garlic for isotropic and pepper for anisotropic."""
    from epr_simfit.model_library import default_components
    from epr_simfit.export import easyspin_script
    from epr_simfit.fitter import fit_spectrum
    from epr_simfit.simulator import simulate_model

    field = np.linspace(345, 358, 400)
    iso = default_components()["nitroxide_14n"].clone()
    y, _ = simulate_model(field, [iso], {"nitroxide_14n": 1.0}, 9.5)
    fit_iso = fit_spectrum(field, y, components=[iso], mw_frequency_GHz=9.5, max_nfev=30)
    assert "garlic(" in easyspin_script(fit_iso, 9.5)

    field2 = np.linspace(280, 350, 500)
    cu = default_components()["cu2_axial"].clone()
    y2, _ = simulate_model(field2, [cu], {"cu2_axial": 1.0}, 9.5, n_orientations=300)
    fit_cu = fit_spectrum(field2, y2, components=[cu], mw_frequency_GHz=9.5, max_nfev=20, n_orientations=300)
    script = easyspin_script(fit_cu, 9.5)
    assert "pepper(" in script
    assert "Sys1.g = [" in script  # g-tensor emitted


def test_fit_reports_uncertainties():
    """Fitted parameter table includes finite weight standard errors."""
    from epr_simfit.model_library import default_components
    from epr_simfit.fitter import fit_spectrum
    from epr_simfit.simulator import simulate_model

    field = np.linspace(345, 358, 800)
    c = default_components()["nitroxide_14n"].clone()
    y, _ = simulate_model(field, [c], {"nitroxide_14n": 1.0}, 9.5)
    yn = y + np.random.default_rng(1).normal(0, 0.01, y.size)
    fit = fit_spectrum(field, yn, components=[c], mw_frequency_GHz=9.5, max_nfev=100)
    assert "std_error" in fit.parameters.columns
    we = fit.weight_errors.get("nitroxide_14n")
    assert we is not None and np.isfinite(we)

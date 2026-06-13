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

"""Physical laws that any correct result must obey."""
import math

import pytest

from cases import CASES, LOSSLESS
from mdrun import MD2D, MD3D, run, h_x, h_xy, sine


@pytest.mark.parametrize("name", LOSSLESS)
def test_energy_conservation(name, tmp_path):
    """Lossless structures: reflected + transmitted efficiencies add up to 1.

    The tolerance reflects the modest NS of the fast cases; the error falls as
    NS grows (see test_energy_improves_with_NS)."""
    program, params, profile, args = CASES[name]
    res = run(program, tmp_path, params, profile, args)
    assert res.returncode == 0
    assert abs(res.sum_eff - 1) < 2e-4


def test_energy_improves_with_NS(tmp_path):
    """TE: the energy balance is a good indicator of the NS convergence. (In TM a small floor,
    about 3e-6 here, comes from the finite-N factorisation and does not decrease with NS.)"""
    program, params, profile, args = CASES["md2D_sine_dielectric_TE_RK4"]
    errors = [abs(run(program, tmp_path / str(ns), dict(params, NS=ns), profile).sum_eff - 1) for ns in (5, 20, 80)]
    assert errors[0] > errors[1] > errors[2]
    assert errors[2] < 1e-9


@pytest.mark.parametrize("pola", ["TE", "TM"])
def test_md2D_symmetric_profile_normal_incidence(pola, tmp_path):
    """A profile symmetric in x, lit at normal incidence, diffracts e(n) = e(-n)."""
    res = run(MD2D, tmp_path, dict(n_sub="1.5 + i0.0", L=1.3, h=0.3, theta_i=0.0, pola=pola, N=8, NS=100,
                                   **{"lambda": 0.6}), h_x(sine(64)))
    for rt in "rt":
        o = res.orders(rt)
        for n in o:
            if n > 0:
                assert o[n] == pytest.approx(o[-n], rel=1e-10, abs=1e-15)


def test_md3D_square_symmetry_swaps_polarisation(tmp_path):
    """Square-symmetric surface at normal incidence: rotating the polarisation by 90 degrees
    exchanges the roles of x and y, so e_psi0(nx, ny) = e_psi90(ny, nx)."""
    n = 32
    prof = h_xy([[0.25 * (math.sin(2 * math.pi * i / n) + math.sin(2 * math.pi * j / n)) for i in range(n)]
                 for j in range(n)])
    p = dict(nu_sub="2.0 + i0.0", Lx=1.0, Ly=1.0, h=0.1, theta_i=0.0, phi_i=0.0, Nx=2, Ny=2, NS=10,
             **{"lambda": 0.83})
    a = run(MD3D, tmp_path / "psi0", dict(p, psi=0.0), prof)
    b = run(MD3D, tmp_path / "psi90", dict(p, psi=90.0), prof)
    for rt in "rt":
        oa, ob = a.orders(rt), b.orders(rt)
        for (nx, ny), e in oa.items():
            assert e == pytest.approx(ob[(ny, nx)], rel=1e-9, abs=1e-14)


def test_reciprocity(tmp_path):
    """Echelette grating: efficiency of the pair (theta_i -> theta_n) equals (-theta_n -> -theta_i).

    Reciprocity holds for the exact solution; at N = 20 the truncated problem satisfies it
    to about 1e-3 relative (4e-5 at N = 40, see thesis case 02)."""
    ech = h_x([1.0 - i / 255 for i in range(256)])
    p = dict(n_sub="1.5 + i0.0", L=1500, h=500, pola="TM", N=20, NS=300, **{"lambda": 632.8})
    a = run(MD2D, tmp_path / "A", dict(p, theta_i=5.0), ech)
    ta = dict(zip(a.arrays["N_eff_r"], a.arrays["theta_r"]))[1.0]        # angle of order +1
    b = run(MD2D, tmp_path / "B", dict(p, theta_i=-ta), ech)
    back = [n for n, t in zip(b.arrays["N_eff_r"], b.arrays["theta_r"]) if abs(t + 5.0) < 1e-4]
    assert back, "no order of B at -5 degrees"
    assert b.eff("r", int(back[0])) == pytest.approx(a.eff("r", 1), rel=2e-3)

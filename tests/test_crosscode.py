"""md2D and md3D must agree where their domains overlap."""
import pytest

from mdrun import MD2D, MD3D, run, h_x, h_xy, sine, multi, lamellar_rows, n_xyz


@pytest.mark.parametrize("pol,psi", [("TE", 0.0), ("TM", 90.0)])
@pytest.mark.parametrize("n_sub", ["1.5 + i0.0", "1.3 + i7.1"])
def test_md3D_Ny0_equals_md2D(pol, psi, n_sub, tmp_path):
    """A 1D grating in md3D (Ny = 0, phi = 0) reproduces md2D in TE (psi = 0) and TM (psi = 90)."""
    common = dict(h=0.3, theta_i=20.0, NS=100, **{"lambda": 0.6})
    a = run(MD2D, tmp_path / "2d", dict(common, n_sub=n_sub, L=1.3, pola=pol, N=8), h_x(sine(128)))
    b = run(MD3D, tmp_path / "3d", dict(common, nu_sub=n_sub, Lx=1.3, Ly=0.1, phi_i=0.0, psi=psi, Nx=8, Ny=0),
            h_xy([sine(128)]))
    oa, ob = a.orders("r"), b.orders("r")
    assert {(n, 0) for n in oa} == set(ob)
    for n, e in oa.items():
        assert ob[(n, 0)] == pytest.approx(e, rel=1e-9, abs=1e-14)


def test_RK4_converges_to_ZINVAR_on_binary_structure_TE(tmp_path):
    """TE, lamellar line whose flat top and bottom lie on slice boundaries: with the one-sided
    rule RK4 converges to the exact z-invariant result at fourth order in NS (first order before)."""
    prof = multi(["2.0 + i0.1"], lamellar_rows(128, 0.4))
    p = dict(n_sub="1.5 + i0.0", L=1.3, h=0.3, theta_i=20.0, pola="TE", N=8, **{"lambda": 0.6})
    exact = run(MD2D, tmp_path / "z", dict(p, calcul_method="Z_INVAR", NS=4), prof).orders("r")

    def worst(ns):
        rk = run(MD2D, tmp_path / str(ns), dict(p, NS=ns), prof).orders("r")
        return max(abs(rk[n] - e) / e for n, e in exact.items())
    coarse, fine = worst(100), worst(400)
    assert fine < 1e-7
    assert coarse / fine > 100          # ~4^4


def test_RK4_converges_to_ZINVAR_on_layered_index_map(tmp_path):
    """Index map made of 4 horizontal rows, row boundaries on slice boundaries: fourth order."""
    nmap = n_xyz([[1.0 if abs(i - 32) < 8 + 4 * k else 1.6 for i in range(64)] for k in range(4)])
    p = dict(n_sub="1.5 + i0.0", L=1.3, h=0.2, theta_i=20.0, pola="TE", N=8, **{"lambda": 0.6})
    exact = run(MD2D, tmp_path / "z", dict(p, calcul_method="Z_INVAR", NS=4), nmap).orders("r")

    def worst(ns):
        rk = run(MD2D, tmp_path / str(ns), dict(p, NS=ns), nmap).orders("r")
        return max(abs(rk[n] - e) / e for n, e in exact.items())
    coarse, fine = worst(40), worst(160)
    assert fine < 1e-7
    assert coarse / fine > 100


def test_ZINVAR_and_IMPROVED_RCWA_converge_together_TM(tmp_path):
    """TM, lamellar line: the exact lamellar formulation (Z_INVAR) and the formulation that uses
    the finite-difference normal of the stepped interfaces (IMPROVED_RCWA, the NS -> inf limit
    of RK4) approach the same value as N grows; the second one more slowly."""
    prof = multi(["2.0 + i0.1"], lamellar_rows(1024, 0.4))
    p = dict(n_sub="1.5 + i0.0", L=1.3, h=0.3, theta_i=20.0, pola="TM", NS=20, **{"lambda": 0.6})

    def gap(N):
        z = run(MD2D, tmp_path / ("z%d" % N), dict(p, calcul_method="Z_INVAR", N=N), prof)
        i = run(MD2D, tmp_path / ("i%d" % N), dict(p, calcul_method="IMPROVED_RCWA", N=N), prof)
        return abs(z.eff("r", -1) - i.eff("r", -1)) / z.eff("r", -1)
    assert gap(32) < gap(8) / 2


def test_ZINVAR_needs_enough_slices_for_large_N(tmp_path):
    """Within one slice evanescent orders grow like exp(2 pi N h_slice / L); with enough slices
    the S-matrix algorithm keeps a thick z-invariant layer stable at large N."""
    prof = multi(["2.0 + i0.1"], lamellar_rows(1024, 0.4))
    p = dict(n_sub="1.5 + i0.0", L=1.3, h=0.3, theta_i=20.0, pola="TM", N=64, calcul_method="Z_INVAR",
             **{"lambda": 0.6})
    res = run(MD2D, tmp_path, dict(p, NS=20), prof)
    assert 0 < res.sum_eff < 1

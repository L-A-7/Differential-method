"""Comparisons with exact analytic solutions (flat interfaces and thin films)."""
import cmath
import math

import pytest

from mdrun import MD2D, MD3D, run, h_x, multi, h_xy


def _cos(n, s):
    c = cmath.sqrt(1 - (s / n) ** 2)
    return -c if c.imag < 0 else c


def _r(na, ca, nb, cb, pol):
    return (na * ca - nb * cb) / (na * ca + nb * cb) if pol == "TE" else (nb * ca - na * cb) / (nb * ca + na * cb)


def fresnel_R(n1, n2, theta_deg, pol):
    s = n1 * math.sin(math.radians(theta_deg))
    return abs(_r(n1, _cos(n1, s), n2, _cos(n2, s), pol)) ** 2


def film_R(n1, nf, n2, d, lam, theta_deg, pol):
    s = n1 * math.sin(math.radians(theta_deg))
    c1, cf, c2 = _cos(n1, s), _cos(nf, s), _cos(n2, s)
    r12, r23 = _r(n1, c1, nf, cf, pol), _r(nf, cf, n2, c2, pol)
    ph = cmath.exp(2j * 2 * math.pi / lam * nf * cf * d)
    return abs((r12 + r23 * ph) / (1 + r12 * r23 * ph)) ** 2


def _fmt(n):
    return "%.15g + i%.15g" % (n.real, n.imag)


SUBSTRATES = [complex(1.5, 0), complex(1.3, 7.1), complex(3.9, 0.02)]


@pytest.mark.parametrize("n2", SUBSTRATES)
@pytest.mark.parametrize("pol", ["TE", "TM"])
@pytest.mark.parametrize("theta", [0.0, 35.0, 70.0])
def test_md2D_flat_interface_fresnel(n2, pol, theta, tmp_path):
    """Flat profile, z-invariant propagation: reflectance equals Fresnel's to rounding."""
    res = run(MD2D, tmp_path, dict(n_sub=_fmt(n2), L=1.0, h=0.2, theta_i=theta, pola=pol, N=3, NS=1,
                                   calcul_method="Z_INVAR", **{"lambda": 0.6}), h_x([0.0] * 32))
    assert res.eff("r", 0) == pytest.approx(fresnel_R(1.0, n2, theta, pol), rel=1e-12, abs=1e-14)


@pytest.mark.parametrize("pol", ["TE", "TM"])
def test_md2D_flat_interface_RK4_converges(pol, tmp_path):
    """RK4 through a homogeneous zone converges to Fresnel as the slices get thinner.

    The interface sits on a slice boundary, where the permittivity is discontinuous,
    so the error falls as NS^-2 rather than NS^-4."""
    ref = fresnel_R(1.0, 1.5, 35.0, pol)
    errs = []
    for ns in (100, 1000):
        res = run(MD2D, tmp_path / str(ns), dict(n_sub="1.5 + i0.0", L=1.0, h=0.2, theta_i=35.0, pola=pol, N=3,
                                                 NS=ns, **{"lambda": 0.6}), h_x([0.0] * 32))
        errs.append(abs(res.eff("r", 0) - ref))
    assert errs[1] < 1e-7
    assert errs[0] / errs[1] > 50          # second order at least


@pytest.mark.parametrize("pol", ["TE", "TM"])
@pytest.mark.parametrize("theta", [0.0, 40.0])
def test_md2D_thin_film_airy(pol, theta, tmp_path):
    """Absorbing film on glass, described as a flat one-layer stack: Airy formula."""
    nf = complex(2.1, 0.05)
    res = run(MD2D, tmp_path, dict(n_sub="1.5 + i0.0", L=1.0, h=0.37, theta_i=theta, pola=pol, N=3, NS=1,
                                   calcul_method="Z_INVAR", **{"lambda": 0.6}),
              multi([_fmt(nf)], [(1.0, 0.0)] * 32))
    assert res.eff("r", 0) == pytest.approx(film_R(1.0, nf, 1.5, 0.37, 0.6, theta, pol), rel=1e-12)


@pytest.mark.parametrize("psi", [0.0, 90.0, 45.0])
@pytest.mark.parametrize("n2", [complex(1.5, 0), complex(1.3, 7.1)])
def test_md3D_flat_interface_any_azimuth(psi, n2, tmp_path):
    """md3D, oblique incidence in an arbitrary plane (phi = 30 deg): psi = 0 is TE, psi = 90 is TM,
    and a mixed polarisation reflects the corresponding mixture."""
    res = run(MD3D, tmp_path, dict(nu_sub=_fmt(n2), Lx=1.0, Ly=1.0, h=0.2, theta_i=35.0, phi_i=30.0, psi=psi,
                                   Nx=1, Ny=1, NS=1, calcul_method="Z_INVAR", **{"lambda": 0.6}),
              h_xy([[0.0] * 16 for _ in range(16)]))
    c, s = math.cos(math.radians(psi)) ** 2, math.sin(math.radians(psi)) ** 2
    ref = c * fresnel_R(1.0, n2, 35.0, "TE") + s * fresnel_R(1.0, n2, 35.0, "TM")
    assert res.eff("r", 0, 0) == pytest.approx(ref, rel=1e-11, abs=1e-14)

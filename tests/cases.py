"""Small, fast configurations covering the features of md2D and md3D.

Used by the regression tests (compared with tests/golden/*.json) and by the
physics tests. Each entry: name -> (program, params, profile_text, extra_args).
Keep every case well under a few seconds.
"""
import math

from mdrun import MD2D, MD3D, h_x, sine, multi, lamellar_rows, n_xyz, h_xy, circle_map

AL = "1.3 + i7.1"          # aluminium at 0.6 um
GLASS = "1.5 + i0.0"

_sine64 = h_x(sine(64))
_sine256 = h_x(sine(256))
_echelette = h_x([1.0 - i / 127 for i in range(128)])
_line = multi(["2.0 + i0.1"], lamellar_rows(128, 0.4))
_two_layers = multi(["1.46 + i0.0", "3.9 + i0.02"],
                    [(1.0, 0.6 + 0.3 * math.cos(2 * math.pi * i / 128), 0.3 * (1 + math.sin(2 * math.pi * i / 128)) / 2)
                     for i in range(128)])
_nmap = n_xyz([[1.0 if abs(i - 32) < 8 + 4 * k else 1.6 for i in range(64)] for k in range(4)])


def _sin3d(n):
    return h_xy([[0.25 * (math.sin(2 * math.pi * i / n) + math.sin(2 * math.pi * j / n)) for i in range(n)]
                 for j in range(n)])


def _pyramid(n):
    return h_xy([[1.0 - max(abs(2 * i / n - 1), abs(2 * j / n - 1)) for i in range(n)] for j in range(n)])


_d2 = dict(n_sub=GLASS, L=1.3, h=0.3, theta_i=20.0, N=8, NS=100, **{"lambda": 0.6})
_m2 = dict(n_sub=AL, L=0.8333333, h=0.3, theta_i=10.0, N=10, NS=300, **{"lambda": 0.6})
_d3 = dict(nu_sub=GLASS, Lx=1.0, Ly=1.0, h=0.1, theta_i=0.0, phi_i=0.0, psi=0.0, Nx=2, Ny=2, NS=10,
           **{"lambda": 0.83})

CASES = {
    # --- md2D, sinusoidal interface
    "md2D_sine_dielectric_TE_RK4": (MD2D, dict(_d2, pola="TE"), _sine64, ()),
    "md2D_sine_dielectric_TM_RK4": (MD2D, dict(_d2, pola="TM"), _sine64, ()),
    "md2D_sine_metal_TE_RK4":      (MD2D, dict(_m2, pola="TE"), _sine256, ()),
    "md2D_sine_metal_TM_RK4":      (MD2D, dict(_m2, pola="TM"), _sine256, ()),
    "md2D_sine_metal_TM_ZINVAR":   (MD2D, dict(_m2, pola="TM", calcul_method="Z_INVAR"), _sine256, ()),
    "md2D_sine_metal_TM_IMPROVED_RCWA": (MD2D, dict(_m2, pola="TM", calcul_method="IMPROVED_RCWA"), _sine256, ()),
    "md2D_sine_normal_TM":         (MD2D, dict(_d2, pola="TM", theta_i=0.0), _sine64, ()),
    "md2D_sine_hamming_TM":        (MD2D, dict(_d2, pola="TM", fft_filter="HAMMING"), _sine64, ()),   # filter acts in TM only
    "md2D_line_smoothing_TM":      (MD2D, dict(_d2, pola="TM", NS=200, smoothing=1, l_smooth=0.02), _line, ()),
    "md2D_sine_auto_TE":           (MD2D, dict(_d2, pola="TE", N="AUTO", NS="AUTO"), _sine64, ()),
    "md2D_echelette_TM":           (MD2D, dict(_d2, pola="TM", L=1.5, h=0.5, theta_i=5.0, N=12, NS=150), _echelette, ()),
    # --- md2D, stacks and index maps
    "md2D_line_TE_ZINVAR":         (MD2D, dict(_d2, pola="TE", calcul_method="Z_INVAR", NS=4), _line, ()),
    "md2D_line_TM_ZINVAR":         (MD2D, dict(_d2, pola="TM", calcul_method="Z_INVAR", NS=4), _line, ()),
    "md2D_line_TM_RK4":            (MD2D, dict(_d2, pola="TM", NS=200), _line, ()),
    "md2D_two_layers_TM_RK4":      (MD2D, dict(_d2, pola="TM", h=0.4, NS=200), _two_layers, ()),
    "md2D_index_map_TE":           (MD2D, dict(_d2, pola="TE", h=0.2, NS=40), _nmap, ()),
    # --- md3D
    "md3D_sin3D_psi0":             (MD3D, dict(_d3), _sin3d(32), ()),
    "md3D_sin3D_psi90":            (MD3D, dict(_d3, psi=90.0), _sin3d(32), ()),
    "md3D_pyramid_oblique":        (MD3D, dict(nu_sub=GLASS, Lx=1.5, Ly=1.0, h=0.25, theta_i=30.0, phi_i=45.0, psi=90.0,
                                                Nx=2, Ny=2, NS=10, **{"lambda": 1.533}), _pyramid(64), ()),
    "md3D_holes_ZINVAR":           (MD3D, dict(nu_sub=GLASS, Lx=1.0001, Ly=1.0001, h=0.05, psi=45.0, Nx=2, Ny=2, NS=1,
                                                calcul_method="Z_INVAR", **{"lambda": 0.5}),
                                    circle_map(64, 0.25, "1.0 + i0.0", "1.75 + i1.5"), ()),
    "md3D_conical_lamellar":       (MD3D, dict(nu_sub=GLASS, Lx=3.1, Ly=0.1, h=1.0, theta_i=30.0, phi_i=45.0, psi=45.0,
                                                Nx=10, Ny=0, NS=40, **{"lambda": 1.0}),
                                    h_xy([[1.0 if i < 128 else 0.0 for i in range(256)]]), ()),
}

# Cases with lossless materials: the efficiencies must add up to 1.
LOSSLESS = [k for k in CASES if "metal" not in k and "line" not in k and "two_layers" not in k
            and "holes" not in k and "auto" not in k]

NEAR_FIELD_CASE = (MD2D, dict(_d2, pola="TM", calcul_type="NEAR_FIELD", verbosity=0, N=6, NS=10), h_x(sine(32)), ())

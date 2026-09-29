"""Validation tier: the thesis tables recomputed at reduced resolution (`make test-validation`, ~2 min).

The reference is always the thesis (validation_cases/*/data.csv). Cost grows as N^3·NS
(md2D) and ((2N+1)^2)^3·NS (md3D), so each case runs at a reduced N or NS chosen to take seconds:

- where the thesis tabulates convergence against N (cases 05, 09/10, 12, 13) the run is compared
  with the thesis value *at the same N*, which is a tight test;
- where only converged values exist (01-04, 06, 11) the tolerance covers the measured
  convergence gap at the reduced N (thresholds are 2-3 x the deviation measured in Sept. 2026).

Dielectrics: 02, 03, 06, 09/10, 11. Metals and absorbing media: 01, 04, 05, 12 (3D), 13.
Full-resolution reproductions of every table: validation_cases/*/run.sh.
"""
import csv
import os
import re
import shutil
import subprocess

import pytest

from mdrun import MD2D, MD3D, REPO, read_arrays, Result

pytestmark = pytest.mark.validation

TVC = os.path.join(REPO, "validation_cases")
MKPROFILE = os.path.join(TVC, "tools", "mkprofile.py")
SIN2048 = "md2D/sin_2048.txt"


def case_dir(prefix):
    return os.path.join(TVC, next(d for d in os.listdir(TVC) if d.startswith(prefix)))


def table(prefix):
    with open(os.path.join(case_dir(prefix), "data.csv")) as f:
        return list(csv.DictReader(f))


def run_case(prefix, workdir, program, args=(), profile=None, generate=None, name="res"):
    """Run a validation case's param.txt in workdir with command-line overrides."""
    workdir = str(workdir)
    os.makedirs(workdir, exist_ok=True)
    shutil.copy(os.path.join(case_dir(prefix), "param.txt"), workdir)
    if profile:
        os.symlink(os.path.join(REPO, profile), os.path.join(workdir, os.path.basename(profile)))
    if generate:
        subprocess.run(["python3", MKPROFILE] + generate, cwd=workdir, check=True)
    proc = subprocess.run([program, "-param", "param.txt", "-verbosity", "0", "-profile_name", name]
                          + [str(a) for a in args], cwd=workdir, capture_output=True, text=True, errors="replace")
    assert proc.returncode == 0, proc.stderr[-500:]
    res = Result(proc.returncode, proc.stdout, proc.stderr)
    res.arrays = read_arrays(os.path.join(workdir, name + ".txt"))
    return res


def pair(label):
    a, b = re.search(r"\((-?\d+),\s*(-?\d+)\)", label).groups()
    return int(a), int(b)


# ---------------------------------------------------------------- md2D

def test_case01_metal_TM_matrix_integration(tmp_path):
    """Aluminium-like sinusoidal grating, TM (thesis N=30, NS=1000). N=20, NS=500: ~4 s."""
    res = run_case("case01", tmp_path, MD2D, ["-N", 20, "-NS", 500], profile=SIN2048)
    for r in table("case01"):
        assert res.eff("r", int(r["order"])) == pytest.approx(float(r["e_matrix_integration"]), abs=1.5e-3)


def test_case02_dielectric_reciprocity(tmp_path):
    """Glass echelette, TM, both reciprocal incidences (thesis N=40). N=20, NS=300: ~5 s."""
    gen = ["echelette", "echelette_256.txt", "--nx", "256"]
    a = run_case("case02", tmp_path / "A", MD2D, ["-N", 20, "-NS", 300, "-theta_i", "5.0"], generate=gen)
    b = run_case("case02", tmp_path / "B", MD2D, ["-N", 20, "-NS", 300, "-theta_i", "-30.60"], generate=gen)
    for r in table("case02"):
        assert a.eff("r", int(r["n_A"])) == pytest.approx(float(r["e_n_A__thetai_5.00deg"]), abs=3e-4)
        assert b.eff("r", int(r["n_B"])) == pytest.approx(float(r["e_n_B__thetai_-30.60deg"]), abs=3e-4)
    assert a.eff("r", 1) == pytest.approx(b.eff("r", 1), rel=2e-4)


@pytest.mark.parametrize("pola", ["TE", "TM"])
def test_case03_dielectric_vs_MFS(pola, tmp_path):
    """Dielectric sinusoidal grating vs fictitious sources (thesis N=100). N=30, NS=300: 4-8 s."""
    res = run_case("case03", tmp_path, MD2D, ["-N", 30, "-NS", 300, "-pola", pola], profile=SIN2048)
    for r in table("case03"):
        if r["polarization"] == pola:
            assert res.eff("r", int(r["order"])) == pytest.approx(float(r["e_MFS"]), abs=1e-5)


@pytest.mark.parametrize("pola", ["TE", "TM"])
def test_case04_metal_vs_MFS(pola, tmp_path):
    """Aluminium sinusoidal grating, thesis DM column (N=40; agrees with MFS to 1.5e-3). N=25, NS=700: 5-11 s."""
    res = run_case("case04", tmp_path, MD2D, ["-N", 25, "-NS", 700, "-pola", pola], profile=SIN2048)
    for r in table("case04"):
        if r["polarization"] == pola:
            n = int(r["order"].split("/")[0])
            assert res.eff("r", n) == pytest.approx(float(r["e_DM"]), abs=3e-3)


@pytest.mark.parametrize("N", [5, 10])
def test_case05_metal_convergence_DM_and_RCWA(N, tmp_path):
    """Aluminium grating, TE, convergence table: both columns compared at the same N. ~1-4 s."""
    rows = {int(r["N"]): r for r in table("case05")}
    z = run_case("case05", tmp_path / "z", MD2D, ["-N", N, "-NS", 1000, "-calcul_method", "Z_INVAR"], profile=SIN2048)
    d = run_case("case05", tmp_path / "d", MD2D, ["-N", N, "-NS", 2000], profile=SIN2048)
    assert z.eff("r", 0) == pytest.approx(float(rows[N]["e0_RCWA"]), rel=1e-5)
    assert d.eff("r", 0) == pytest.approx(float(rows[N]["e0_differential_method"]), rel=1.5e-4)


# ---------------------------------------------------------------- md3D

def test_case06_dielectric_conical_vs_gsolver(tmp_path):
    """Lamellar glass grating, conical incidence (thesis Nx=50). Nx=20, NS=200: ~10 s.
    Within 3.2e-3 of GSolver since the one-sided slice-boundary rule (1.3e-2 before)."""
    res = run_case("case06", tmp_path, MD3D, ["-Nx", 20, "-NS", 200], profile="md3D/carre_256_1D.txt")
    for r in table("case06"):
        rt = "r" if r["direction"] == "reflection" else "t"
        assert res.eff(rt, -int(r["order"]), 0) == pytest.approx(float(r["e_GSolver"]), rel=1e-2)


@pytest.mark.parametrize("N", [2, 3])
def test_case09_case10_dielectric_pyramid(N, tmp_path):
    """Glass pyramid, convergence tables (reflection and transmission) at the same N. 1.5-3 s."""
    res = run_case("case08", tmp_path, MD3D, ["-Nx", N, "-Ny", N],
                   generate=["pyramid", "pyramid_512x512.txt", "--n", "512"])
    for prefix, rt in (("case09", "r"), ("case10", "t")):
        for r in table(prefix):
            if r["method"] == "differential method" and int(r["N"]) == N:
                for col in (c for c in r if c.startswith(rt + "(")):
                    a, b = pair(col)
                    assert res.eff(rt, b, a) == pytest.approx(float(r[col]), rel=1.5e-3), col


def test_case11_dielectric_biperiodic_sine(tmp_path):
    """Biperiodic sinusoid vs variation of boundaries (thesis N=15). N=3: ~5 s."""
    res = run_case("case11", tmp_path, MD3D, ["-Nx", 3, "-Ny", 3], profile="md3D/pr_sin3D_256x256.txt")
    for r in table("case11"):
        assert res.eff("r", *pair(r["quantity"])) == pytest.approx(float(r["differential_method"]), rel=2e-4)


@pytest.mark.parametrize("N", [0, 1, 2])
def test_case12_metal_deep_aluminium(N, tmp_path):
    """Deep aluminium biperiodic grating, convergence table at the same N (NS=250). 4-7 s."""
    rows = {int(r["N"]): r for r in table("case12")}
    res = run_case("case12", tmp_path, MD3D, ["-Nx", N, "-Ny", N], profile="md3D/pr_sin3D_256x256.txt")
    for col in ("e(-1,0)", "e(0,-1)", "e(0,0)"):
        if rows[N][col]:
            assert res.eff("r", *pair(col)) == pytest.approx(float(rows[N][col]), rel=5e-4), col


@pytest.mark.parametrize("N", [0, 1, 2, 3, 4])
def test_case13_absorbing_holes(N, tmp_path):
    """Holes in an absorbing layer vs Schuster et al., convergence table at the same N. < 1 s each."""
    rows = {int(r["N"]): r for r in table("case13")}
    res = run_case("case13", tmp_path, MD3D, ["-Nx", N, "-Ny", N, "-profile_file", "circle_R025_512x512.txt"],
                   generate=["circle", "circle_R025_512x512.txt", "--n", "512", "--R", "0.25", "--L", "1.0"])
    assert res.eff("r", 0, 0) == pytest.approx(float(rows[N]["md3D_NV_radial"]), abs=1.5e-4)

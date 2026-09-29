#!/usr/bin/env python3
"""Generate the profile files needed by the thesis validation cases.

The original profile files of several thesis cases were lost; the shapes are
simple enough to be regenerated exactly from their analytic description.
Profiles are normalised by md2D/md3D to the height h given in the parameter
file, so only the shape matters here.

Usage:
  mkprofile.py echelette OUT [--nx 256]
      md2D H_X blazed (single-slope) profile, one vertical wall per period.
      Slope oriented so that md2D reproduces thesis Table 2.2 (the old md1D
      file echelle.txt had the opposite orientation).
  mkprofile.py pyramid OUT [--n 512]
      md3D H_XY rectangular-base pyramid filling the whole period,
      z = 1 - max(|2x/Lx - 1|, |2y/Ly - 1|), x = j*Lx/n (apex on a sample).
  mkprofile.py circle OUT [--n 1024] [--R 0.25] [--L 1.0] [--n1 1.0+i0.0] [--n2 1.75+i1.5]
      md3D N_XY_ZINVAR index map: disk of radius R (index n1) in a square
      cell of side L (index n2). Same sampling as md3D/utils/profilGen2D.c
      CIRCLE, whose centre also matches md3D's hard-coded radial normal field.
"""
import argparse


def echelette(out, nx):
    with open(out, "w") as f:
        f.write("type_profil = H_X\nN_x = %d\nprofil =\n" % nx)
        for i in range(nx):
            f.write("%.15g\n" % (1.0 - i / (nx - 1)))


def pyramid(out, n):
    with open(out, "w") as f:
        f.write("profile_type = H_XY\nNprx = %d\nNpry = %d\nprofil =\n" % (n, n))
        for iy in range(n):
            y = iy / n
            f.write(" ".join("%.9e" % (1.0 - max(abs(2 * ix / n - 1), abs(2 * y - 1)))
                             for ix in range(n)) + "\n")


def circle(out, n, R, L, n1, n2):
    c = (n - 1) / 2
    with open(out, "w") as f:
        f.write("n1 = %s\nn2 = %s\n" % (n1, n2))
        f.write("profile_type = N_XY_ZINVAR\nNprx = %d\nNpry = %d\nNprz = 1\nn_xyz = \n" % (n, n))
        for iy in range(n):
            dy2 = ((iy - c) * L / n) ** 2
            f.write(" ".join("1" if ((ix - c) * L / n) ** 2 + dy2 <= R * R else "2"
                             for ix in range(n)) + "\n")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("kind", choices=["echelette", "pyramid", "circle"])
    ap.add_argument("out")
    ap.add_argument("--nx", type=int, default=256)
    ap.add_argument("--n", type=int, default=None)
    ap.add_argument("--R", type=float, default=0.25)
    ap.add_argument("--L", type=float, default=1.0)
    ap.add_argument("--n1", default="1.0 + i0.0")
    ap.add_argument("--n2", default="1.75 + i1.5")
    a = ap.parse_args()
    if a.kind == "echelette":
        echelette(a.out, a.nx)
    elif a.kind == "pyramid":
        pyramid(a.out, a.n or 512)
    else:
        circle(a.out, a.n or 1024, a.R, a.L, a.n1, a.n2)

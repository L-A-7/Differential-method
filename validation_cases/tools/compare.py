#!/usr/bin/env python3
"""Compare md2D/md3D results produced by a case's run.sh with the thesis values.

Usage: compare.py CASE_DIR        (reads CASE_DIR/data.csv and CASE_DIR/run/*.txt)

Order-label conventions (checked against the thesis tables):
- md2D: thesis order n == md2D order n (cases 01-05); case 06 (md3D, Ny=0)
  thesis order n == md3D nx = -n.
- case 08, 11, 12: thesis (a,b) == md3D (nx=a, ny=b).
- case 09, 10:     thesis (a,b) == md3D (nx=b, ny=a).
"""
import csv
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mdres import eff, read_arrays  # noqa: E402

CASE = os.path.abspath(sys.argv[1])
RUN = os.path.join(CASE, "run")
NAME = os.path.basename(CASE)[:6]
ROWS = list(csv.DictReader(open(os.path.join(CASE, "data.csv"))))


def res(fname):
    p = os.path.join(RUN, fname)
    return p if os.path.exists(p) else None


def get(fname, rt, *order):
    p = res(fname)
    return eff(p, rt, *order) if p else None


def line(label, now, ref, ref_name="thesis"):
    if now is None:
        print("  %-28s %14s   %s=%s" % (label, "(not run)", ref_name, ref))
        return
    try:
        r = float(ref)
        d = "abs=%+.2e rel=%+.2e" % (now - r, (now / r - 1) if r else float("nan"))
    except (TypeError, ValueError):
        d = ""
    print("  %-28s %14.8g   %s=%-12s %s" % (label, now, ref_name, ref, d))


def nfiles(prefix):
    """N values available for files named PREFIX<N>.txt in run/."""
    out = []
    for f in os.listdir(RUN) if os.path.isdir(RUN) else []:
        m = re.match(re.escape(prefix) + r"(\d+)\.txt$", f)
        if m:
            out.append(int(m.group(1)))
    return sorted(out)


def pair(label):
    a, b = re.search(r"\((-?\d+),\s*(-?\d+)\)", label).groups()
    return int(a), int(b)


print("== %s  (reproduced vs thesis)" % os.path.basename(CASE))

if NAME == "case01":
    for r in ROWS:
        line("r order %s" % r["order"], get("c01.txt", "r", int(r["order"])), r["e_matrix_integration"])

elif NAME == "case02":
    for r in ROWS:
        line("A (theta_i=5.00) n=%s" % r["n_A"], get("c02_A.txt", "r", int(r["n_A"])), r["e_n_A__thetai_5.00deg"])
    for r in ROWS:
        line("B (theta_i=-30.60) n=%s" % r["n_B"], get("c02_B.txt", "r", int(r["n_B"])), r["e_n_B__thetai_-30.60deg"])

elif NAME == "case03":
    for r in ROWS:
        now = get("c03_%s.txt" % r["polarization"], "r", int(r["order"]))
        line("%s r order %s" % (r["polarization"], r["order"]), now, r["e_DM"], "DM")
        line("", now, r["e_MFS"], "MFS")

elif NAME == "case04":
    for r in ROWS:
        n = int(r["order"].split("/")[0])
        now = get("c04_%s.txt" % r["polarization"], "r", n)
        line("%s r order %s" % (r["polarization"], r["order"]), now, r["e_DM"], "DM")
        line("", now, r["e_MFS"], "MFS")

elif NAME == "case05":
    for r in ROWS:
        N = int(r["N"])
        line("N=%d DM (RK4)" % N, get("c05_RK4_N%d.txt" % N, "r", 0), r["e0_differential_method"], "DM")
        line("N=%d RCWA (Z_INVAR)" % N, get("c05_ZINVAR_N%d.txt" % N, "r", 0), r["e0_RCWA"], "RCWA")

elif NAME == "case06":
    for r in ROWS:
        rt = "r" if r["direction"] == "reflection" else "t"
        now = get("c06.txt", rt, -int(r["order"]), 0)
        line("%s order %s" % (r["direction"], r["order"]), now, r["e_DM"], "DM")
        line("", now, r["e_GSolver"], "GSolver")

elif NAME == "case07":
    print("  (absolute times depend on the machine; compare the growth ~N^6)")
    t1 = None
    for r in ROWS:
        N = int(r["N"])
        p = res("c07_N%d.txt" % N)
        now = read_arrays(p)["Calcul_duration"][0] if p else None
        if N == 1:
            t1 = (now, float(r["time_seconds"]))
        line("N=%d time [s]" % N, now, r["time_seconds"])
        if now and t1 and t1[0]:
            print("  %-28s %14.4g   thesis=%.4g" % ("   t(N)/t(1)", now / t1[0], float(r["time_seconds"]) / t1[1]))

elif NAME == "case08":
    Ns = nfiles("c08_N")
    if not Ns:
        print("  (not run)")
    else:
        f = "c08_N%d.txt" % Ns[-1]
        print("  using N=%d (thesis: N=15)" % Ns[-1])
        for r in ROWS:
            rt = r["quantity"][2]
            line(r["quantity"], get(f, rt, *pair(r["quantity"])), r["md3D"], "md3D")
            line("", get(f, rt, *pair(r["quantity"])), r["Granet_Chandezon_1998"], "Chandezon")

elif NAME in ("case09", "case10"):
    run8 = os.path.join(os.path.dirname(CASE), "case08_3d_pyramid_dielectric_multi_method", "run")
    RUN = run8
    rt = "r" if NAME == "case09" else "t"
    cols = [c for c in ROWS[0] if c.startswith(rt + "(")]
    for r in ROWS:
        if r["method"] != "differential method":
            continue
        N = int(r["N"])
        for c in cols:
            a, b = pair(c)
            line("N=%d %s" % (N, c), get("c08_N%d.txt" % N, rt, b, a), r[c])

elif NAME in ("case11", "case12"):
    cols = [c for c in ROWS[0] if c.startswith("e(")]
    if NAME == "case11":
        Ns = nfiles("c11_N")
        f = "c11_N%d.txt" % Ns[-1] if Ns else "none"
        print("  using N=%s (thesis: N=15)" % (Ns[-1] if Ns else "-"))
        for r in ROWS:
            line(r["quantity"], get(f, "r", *pair(r["quantity"])), r["differential_method"], "DM")
            line("", get(f, "r", *pair(r["quantity"])), r["variation_of_boundaries"], "MVB")
    else:
        for r in ROWS:
            N = int(r["N"])
            for c in cols:
                if r[c]:
                    line("N=%d %s" % (N, c), get("c12_N%d.txt" % N, "r", *pair(c)), r[c])

elif NAME == "case13":
    for r in ROWS:
        N = int(r["N"])
        now = get("c13_N%d.txt" % N, "r", 0, 0)
        line("N=%d r(0,0)" % N, now, r["md3D_NV_radial"], "md3D")
        line("", now, r["Schuster_NV_radial"], "Schuster")

else:
    sys.exit("unknown case %s" % CASE)

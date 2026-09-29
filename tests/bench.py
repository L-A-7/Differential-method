#!/usr/bin/env python3
"""Timing benchmark for the optimisation work (`make bench`).

Runs a fixed set of configurations, prints the wall-clock and CPU time of each, and appends
one line per run to tests/bench_history.jsonl (commit, date, host, load, times) so that
successive builds can be compared. CPU time (user + system of the solver process) is the
more robust figure on a shared machine; wall time is what users see. For a clean reference,
run on an idle machine (the load before the run is recorded).

    python3 tests/bench.py [--quick] [--repeat N]
"""
import argparse
import datetime
import json
import math
import os
import platform
import resource
import subprocess
import sys
import tempfile
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mdrun import MD2D, MD3D, REPO, run, h_x, sine, h_xy  # noqa: E402


def _pyramid(n):
    return h_xy([[1.0 - max(abs(2 * i / n - 1), abs(2 * j / n - 1)) for i in range(n)] for j in range(n)])


_sine2048 = h_x(sine(2048))
_pyr = _pyramid(256)
_al = dict(n_sub="1.3 + i7.1", L=0.8333333, h=0.8, theta_i=0.0, **{"lambda": 0.6})
_pyr_p = dict(nu_sub="1.5 + i0.0", Lx=1.5, Ly=1.0, h=0.25, theta_i=30.0, phi_i=45.0, psi=90.0, NS=20,
              **{"lambda": 1.533})

BENCHMARKS = [
    # name, program, params, profile, in quick set
    ("md2D TE  N=40  NS=1000 (Al)", MD2D, dict(_al, pola="TE", N=40, NS=1000), _sine2048, True),
    ("md2D TM  N=40  NS=1000 (Al)", MD2D, dict(_al, pola="TM", N=40, NS=1000), _sine2048, True),
    ("md2D TE  N=100 NS=300 (glass)", MD2D, dict(n_sub="1.732 + i0.0", L=3.9, h=0.5, theta_i=30.0, pola="TE", N=100,
                                                 NS=300, **{"lambda": 1.0}), _sine2048, False),
    ("md2D TM  N=40  Z_INVAR NS=200", MD2D, dict(_al, pola="TM", N=40, NS=200, calcul_method="Z_INVAR"), _sine2048, True),
    ("md3D pyramid N=4 NS=20", MD3D, dict(_pyr_p, Nx=4, Ny=4), _pyr, True),
    ("md3D pyramid N=6 NS=20", MD3D, dict(_pyr_p, Nx=6, Ny=6), _pyr, False),
]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--quick", action="store_true", help="only the short benchmarks (about 2 min)")
    ap.add_argument("--repeat", type=int, default=1, help="runs per benchmark; the minimum time is kept")
    args = ap.parse_args()

    commit = subprocess.run(["git", "describe", "--always", "--dirty"], cwd=REPO,
                            capture_output=True, text=True).stdout.strip()
    times, cpu_times = {}, {}
    load_before = os.getloadavg()[0]
    print("%-34s %9s %9s   %s" % ("benchmark", "wall [s]", "cpu [s]", "sum_eff"))
    with tempfile.TemporaryDirectory() as tmp:
        for k, (name, prog, params, prof, quick) in enumerate(BENCHMARKS):
            if args.quick and not quick:
                continue
            best, best_cpu, res = math.inf, math.inf, None
            for rep in range(args.repeat):
                r0 = resource.getrusage(resource.RUSAGE_CHILDREN)
                t0 = time.perf_counter()
                res = run(prog, os.path.join(tmp, "%d_%d" % (k, rep)), params, prof)
                best = min(best, time.perf_counter() - t0)
                r1 = resource.getrusage(resource.RUSAGE_CHILDREN)
                best_cpu = min(best_cpu, (r1.ru_utime - r0.ru_utime) + (r1.ru_stime - r0.ru_stime))
            if res.returncode != 0:
                print("%-34s FAILED: %s" % (name, res.stderr[-200:]))
                continue
            times[name] = round(best, 3)
            cpu_times[name] = round(best_cpu, 3)
            print("%-34s %9.2f %9.2f   %.10f" % (name, best, best_cpu, res.sum_eff))
    entry = dict(commit=commit, date=datetime.datetime.now().isoformat(timespec="seconds"),
                 host=platform.node(), cpu=platform.processor() or platform.machine(),
                 load=round(load_before, 2), times=times, cpu_times=cpu_times)
    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "bench_history.jsonl"), "a") as f:
        f.write(json.dumps(entry) + "\n")
    if load_before > 0.5:
        print("\nNote: load average %.1f before the run - wall times are not those of an idle machine;"
              " compare CPU times" % load_before)


if __name__ == "__main__":
    main()

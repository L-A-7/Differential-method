#!/usr/bin/env python3
"""Parse md2D / md3D result files.

md3D:  eff(path, 'r', nx, ny)   md2D:  eff(path, 'r', n)
CLI:   mdres.py file r|t nx [ny]
"""
import re
import sys


def read_arrays(path):
    """Return dict name -> list of floats for every 'name = v v v ...' block."""
    arrays, cur = {}, None
    for line in open(path, encoding="latin1"):
        m = re.match(r"^\s*([A-Za-z_]\w*)\s*[=:]\s*(.*)$", line)
        if m:
            cur = m.group(1)
            arrays[cur] = []
            rest = m.group(2)
        elif cur is not None:
            rest = line
        else:
            continue
        for tok in rest.split():
            try:
                arrays[cur].append(float(tok))
            except ValueError:
                cur = None  # scalar text field like 'STD', stop collecting
                break
    return arrays


def eff(path, rt, nx, ny=None):
    a = read_arrays(path)
    if ny is None:  # md2D
        for n, e in zip(a["N_eff_" + rt], a["eff_" + rt]):
            if n == nx:
                return e
    else:
        for x, y, e in zip(a["nx_eff_" + rt], a["ny_eff_" + rt], a["eff_" + rt]):
            if x == nx and y == ny:
                return e
    return None


if __name__ == "__main__":
    args = sys.argv[1:]
    print(eff(args[0], args[1], *[int(v) for v in args[2:]]))

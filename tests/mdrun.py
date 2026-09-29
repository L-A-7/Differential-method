"""Helpers to run md2D / md3D on generated inputs and read their results.

Each run happens in its own directory, because the programs write
<profile_name>.txt into the current directory.
"""
import math
import os
import re
import subprocess
from dataclasses import dataclass, field

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MD2D = os.environ.get("MD2D", os.path.join(REPO, "md2D", "md2D"))
MD3D = os.environ.get("MD3D", os.path.join(REPO, "md3D", "md3D"))


# ---------------------------------------------------------------- parameter files

MD2D_DEFAULTS = dict(
    n_super="1.0 + i0.0", i_field_mode="PLANE_WAVE", calcul_type="STD",
    calcul_method="RK4", fft_filter="NONE", smoothing=0, verbosity=1,
    delta_h=0.001, pola="TE", theta_i=0.0,
)
MD3D_DEFAULTS = dict(
    nu_super="1.0 + i0.0", i_field_mode="PLANE_WAVE", calcul_type="STD",
    calcul_method="RK4", verbosity=1, delta_h=0.001,
    theta_i=0.0, phi_i=0.0, psi=0.0,
)


def param_text(defaults, params):
    p = dict(defaults)
    p.update(params)
    return "".join("%s = %s\n" % (k, v) for k, v in p.items())


# ---------------------------------------------------------------- structure files

def h_x(values):
    """md2D single interface z = h(x)."""
    return ("type_profil = H_X\nN_x = %d\nprofil =\n" % len(values)
            + "".join("%.15g\n" % v for v in values))


def sine(n, phase=0.0):
    return [0.5 + 0.5 * math.sin(2 * math.pi * i / n + phase) for i in range(n)]


def multi(indices, rows):
    """md2D stack: indices = ['1.5 + i0.0', ...], rows = [(z1, z2, ...), ...] top to bottom."""
    head = "".join("n%d = %s\n" % (k + 1, v) for k, v in enumerate(indices))
    head += "type_profil = MULTICOUCHES\nN_x = %d\nN_layers = %d\nprofil =\n" % (len(rows), len(indices))
    return head + "".join(" ".join("%.15g" % z for z in r) + "\n" for r in rows)


def lamellar_rows(n, fill, h=1.0):
    """Rows of a one-layer stack: a line of layer material over a fraction `fill` of the period."""
    return [(h, 0.0) if i < round(fill * n) else (h, h) for i in range(n)]


def n_xyz(re_rows, im_rows=None):
    """md2D index map, rows from top (z = h) to bottom."""
    nz, nx = len(re_rows), len(re_rows[0])
    s = "type_profil = N_XYZ\nN_x = %d\nN_z = %d\nRe_n_xyz =\n" % (nx, nz)
    s += "".join(" ".join("%.15g" % v for v in r) + "\n" for r in re_rows)
    if im_rows is not None:
        s += "Im_n_xyz =\n" + "".join(" ".join("%.15g" % v for v in r) + "\n" for r in im_rows)
    return s


def h_xy(z):
    """md3D single interface, z[ny][nx]."""
    s = "profile_type = H_XY\nNprx = %d\nNpry = %d\nprofil =\n" % (len(z[0]), len(z))
    return s + "".join(" ".join("%.15g" % v for v in row) + "\n" for row in z)


def circle_map(n, radius, n_in, n_out, period=1.0):
    c = (n - 1) / 2
    s = "n1 = %s\nn2 = %s\nprofile_type = N_XY_ZINVAR\nNprx = %d\nNpry = %d\nNprz = 1\nn_xyz =\n" % (n_in, n_out, n, n)
    for j in range(n):
        s += " ".join("1" if ((i - c) * period / n) ** 2 + ((j - c) * period / n) ** 2 <= radius ** 2 else "2"
                      for i in range(n)) + "\n"
    return s


# ---------------------------------------------------------------- results

def read_arrays(path):
    """name -> list of floats for every 'name = v v v' block of a result file."""
    arrays, cur = {}, None
    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            m = re.match(r"^\s*([A-Za-z_]\w*)\s*[=:]\s*(.*)$", line)
            if m:
                cur, rest = m.group(1), m.group(2)
                arrays[cur] = []
            elif cur is not None:
                rest = line
            else:
                continue
            for tok in rest.split():
                try:
                    arrays[cur].append(float(tok))
                except ValueError:
                    cur = None
                    break
    return arrays


@dataclass
class Result:
    returncode: int
    stdout: str
    stderr: str
    arrays: dict = field(default_factory=dict)

    def eff(self, rt, nx, ny=None):
        """Efficiency of reflected ('r') or transmitted ('t') order; None if absent."""
        a = self.arrays
        if ny is None and ("N_eff_" + rt) in a:
            for n, e in zip(a["N_eff_" + rt], a["eff_" + rt]):
                if n == nx:
                    return e
            return None
        if ny is None:
            ny = 0
        for x, y, e in zip(a["nx_eff_" + rt], a["ny_eff_" + rt], a["eff_" + rt]):
            if x == nx and y == ny:
                return e
        return None

    def orders(self, rt):
        """{order: efficiency} for propagating orders (md2D: n, md3D: (nx, ny))."""
        a = self.arrays
        if ("N_eff_" + rt) in a:
            return {int(n): e for n, e in zip(a["N_eff_" + rt], a["eff_" + rt])}
        return {(int(x), int(y)): e for x, y, e in zip(a["nx_eff_" + rt], a["ny_eff_" + rt], a["eff_" + rt])
                if math.isfinite(e)}

    @property
    def sum_eff(self):
        return self.arrays["sum_eff"][0]


def run(program, workdir, params, profile, args=(), name="res", defaults=None):
    """Write param.txt and profile.txt in workdir, run program, parse <name>.txt."""
    workdir = str(workdir)
    os.makedirs(workdir, exist_ok=True)
    if defaults is None:
        defaults = MD2D_DEFAULTS if "md2D" in os.path.basename(program) else MD3D_DEFAULTS
    p = dict(params)
    p.setdefault("profile_name", name)
    p.setdefault("profile_file", "profile.txt")
    with open(os.path.join(workdir, "param.txt"), "w") as f:
        f.write(param_text(defaults, p))
    if profile is not None:
        with open(os.path.join(workdir, p["profile_file"]), "w") as f:
            f.write(profile)
    proc = subprocess.run([program, "-param", "param.txt"] + [str(a) for a in args],
                          cwd=workdir, capture_output=True, text=True, errors="replace", timeout=600)
    res = Result(proc.returncode, proc.stdout, proc.stderr)
    out = os.path.join(workdir, p["profile_name"] + ".txt")
    if proc.returncode == 0 and os.path.exists(out):
        res.arrays = read_arrays(out)
    return res


def stdout_arrays(text):
    """Parse labelled arrays printed on stdout (near-field maps)."""
    path = None
    arrays, cur = {}, None
    for line in text.splitlines():
        m = re.match(r"^\s*([A-Za-z_]\w*)\s*=\s*(.*)$", line)
        if m:
            cur, rest = m.group(1), m.group(2)
            arrays[cur] = []
        elif cur is not None:
            rest = line
        else:
            continue
        for tok in rest.split():
            try:
                arrays[cur].append(float(tok))
            except ValueError:
                cur = None
                break
    return arrays

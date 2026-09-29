"""Reading of parameter files, command-line overrides and error handling."""
import pytest

from mdrun import MD2D, MD3D, run, h_x, h_xy, sine

BASE = dict(n_sub="1.5 + i0.0", L=1.3, h=0.3, theta_i=20.0, pola="TE", N=4, NS=20, **{"lambda": 0.6})
PROF = h_x(sine(64))


@pytest.mark.parametrize("text,expected", [
    ("1.3 + i7.6", complex(1.3, 7.6)), ("1.3+i7.6", complex(1.3, 7.6)),
    ("1.3 - i7.6", complex(1.3, -7.6)), ("1.3-i7.6", complex(1.3, -7.6)),
    ("1.3 + i-7.6", complex(1.3, -7.6)), ("1.3 + i 7.6", complex(1.3, 7.6)),
    ("1.3 + 7.6i", complex(1.3, 7.6)), ("1.3 - 7.6*i", complex(1.3, -7.6)),
    ("1.5", complex(1.5, 0.0)), ("1.3e0 + i7.6e0", complex(1.3, 7.6)),
])
def test_complex_syntax(text, expected, tmp_path):
    res = run(MD2D, tmp_path, dict(BASE, n_sub=text), PROF)
    assert res.returncode == 0
    head = next(l for l in open(tmp_path / "res.txt") if l.startswith("n_sub"))
    sign = "-" if expected.imag < 0 else "+"
    assert head.split("=", 1)[1].strip() == "%f %s i%f" % (expected.real, sign, abs(expected.imag))


@pytest.mark.parametrize("text", ["1.3 + 7.6", "abc", ""])
def test_invalid_complex_is_rejected(text, tmp_path):
    res = run(MD2D, tmp_path, dict(BASE, n_sub=text), PROF)
    assert res.returncode != 0
    assert 'can\'t read "n_sub"' in res.stderr


def test_command_line_overrides_file(tmp_path):
    res = run(MD2D, tmp_path, BASE, PROF, args=["-theta_i", "33", "-N", "5"])
    assert res.arrays["theta_i"][0] == pytest.approx(33.0)
    assert res.arrays["N"][0] == 5


def test_command_line_complex(tmp_path):
    res = run(MD2D, tmp_path, BASE, PROF, args=["-sigma0_normed", "0.2 - i0.0"])
    assert res.returncode == 0
    assert "sigma0_normed = 0.2000000000" in res.stdout


def test_missing_smoothing_defaults_to_zero(tmp_path):
    from mdrun import MD2D_DEFAULTS
    defaults = {k: v for k, v in MD2D_DEFAULTS.items() if k != "smoothing"}
    res = run(MD2D, tmp_path, BASE, PROF, defaults=defaults)
    assert res.returncode == 0


@pytest.mark.parametrize("drop", ["L", "lambda", "pola", "fft_filter", "profile_file"])
def test_missing_required_key_fails_cleanly(drop, tmp_path):
    from mdrun import MD2D_DEFAULTS
    params = {k: v for k, v in BASE.items() if k != drop}
    defaults = {k: v for k, v in MD2D_DEFAULTS.items() if k != drop}
    profile = PROF
    if drop == "profile_file":
        params["profile_file"] = "does_not_exist.txt"
        profile = None
    res = run(MD2D, tmp_path, params, profile, defaults=defaults)
    assert res.returncode != 0
    assert "ERROR" in res.stderr


@pytest.mark.parametrize("key,value", [("pola", "XX"), ("calcul_method", "FOO"), ("calcul_type", "BAR")])
def test_invalid_choice_fails_cleanly(key, value, tmp_path):
    res = run(MD2D, tmp_path, dict(BASE, **{key: value}), PROF)
    assert res.returncode != 0
    assert "ERROR" in res.stderr


def test_md3D_missing_key_fails_cleanly(tmp_path):
    res = run(MD3D, tmp_path, dict(nu_sub="1.5", Lx=1.0, h=0.2, Nx=1, Ny=1, NS=2, **{"lambda": 0.6}),
              h_xy([[0.0] * 8] * 8))
    assert res.returncode != 0


def test_profile_name_with_subfolder(tmp_path):
    (tmp_path / "sub").mkdir()
    res = run(MD2D, tmp_path, dict(BASE, profile_name="sub/x"), PROF, name="sub/x")
    assert res.returncode == 0 and (tmp_path / "sub" / "x.txt").exists()

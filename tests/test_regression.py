"""Regression tests: every case of cases.py must reproduce its stored results.

The golden files in tests/golden/ were produced by the code at the commit that
created them (see "commit" inside each file). Regenerate them only after a
change that is *meant* to alter results:  make update-golden
"""
import json
import math
import os
import subprocess

import pytest

from cases import CASES, NEAR_FIELD_CASE
from conftest import ATOL, RTOL
from mdrun import REPO, run, stdout_arrays

GOLDEN = os.path.join(os.path.dirname(os.path.abspath(__file__)), "golden")


def _key(k):
    return "%d,%d" % k if isinstance(k, tuple) else str(k)


def _summary(res):
    return {"r": {_key(k): v for k, v in sorted(res.orders("r").items())},
            "t": {_key(k): v for k, v in sorted(res.orders("t").items())},
            "sum_eff": res.sum_eff}


def _commit():
    try:
        return subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=REPO,
                              capture_output=True, text=True).stdout.strip()
    except OSError:
        return "unknown"


def _close(a, b):
    return abs(a - b) <= max(ATOL, RTOL * max(abs(a), abs(b)))


def _check_or_update(name, data, update):
    path = os.path.join(GOLDEN, name + ".json")
    if update:
        os.makedirs(GOLDEN, exist_ok=True)
        with open(path, "w") as f:
            json.dump(dict(data, commit=_commit()), f, indent=1, sort_keys=True)
        pytest.skip("golden file written")
    if not os.path.exists(path):
        pytest.fail("no golden file %s (run `make update-golden`)" % path)
    with open(path) as f:
        return json.load(f)


@pytest.mark.parametrize("name", sorted(CASES))
def test_case_matches_golden(name, tmp_path, update_golden):
    program, params, profile, args = CASES[name]
    res = run(program, tmp_path, params, profile, args)
    assert res.returncode == 0, res.stderr[-500:]
    got = _summary(res)
    ref = _check_or_update(name, got, update_golden)
    for rt in ("r", "t"):
        assert got[rt].keys() == ref[rt].keys(), "%s: different set of propagating orders" % rt
        bad = {k: (got[rt][k], ref[rt][k]) for k in ref[rt] if not _close(got[rt][k], ref[rt][k])}
        assert not bad, "%s orders differ (now, golden): %s" % (rt, bad)
    assert _close(got["sum_eff"], ref["sum_eff"])


def test_near_field_matches_golden(tmp_path, update_golden):
    program, params, profile, args = NEAR_FIELD_CASE
    res = run(program, tmp_path, params, profile, args)
    assert res.returncode == 0, res.stderr[-500:]
    maps = {k: v for k, v in stdout_arrays(res.stdout).items() if k[:2] in ("Re", "Im")}
    assert set(maps) == {"ReEx", "ImEx", "ReEz", "ImEz", "ReHpy", "ImHpy"}
    got = dict(_summary(res), maps=maps)
    ref = _check_or_update("md2D_near_field_TM", got, update_golden)
    for k, values in ref["maps"].items():
        assert len(got["maps"][k]) == len(values)
        scale = max(abs(v) for v in values)
        worst = max(abs(a - b) for a, b in zip(got["maps"][k], values))
        assert worst <= max(1e-12, RTOL * scale), "%s differs by %.3g" % (k, worst)
        assert all(math.isfinite(v) for v in got["maps"][k])

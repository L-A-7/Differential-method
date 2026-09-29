"""Builds and runs the C unit tests of the shared library code (md_libs)."""
import os
import subprocess

from mdrun import REPO

HERE = os.path.dirname(os.path.abspath(__file__))


def test_md_io_utils(tmp_path):
    exe = tmp_path / "test_md_io_utils"
    build = subprocess.run(
        ["gcc", "-std=gnu99", "-O2", "-I", os.path.join(REPO, "md2D"), "-I", os.path.join(REPO, "md_libs"),
         "-o", str(exe), os.path.join(HERE, "c", "test_md_io_utils.c"),
         os.path.join(REPO, "md_libs", "md_io_utils.c"), "-lm"],
        capture_output=True, text=True)
    assert build.returncode == 0, build.stderr
    res = subprocess.run([str(exe)], cwd=tmp_path, capture_output=True, text=True)
    assert res.returncode == 0, res.stdout + res.stderr

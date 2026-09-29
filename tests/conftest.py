import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mdrun import MD2D, MD3D  # noqa: E402

# Relative tolerance of the regression tests. Results are bit-identical between
# -O0 and -O3 builds on x86-64 today; 1e-9 leaves room for a BLAS change or a
# reordering of sums while catching any real change of the physics.
RTOL = float(os.environ.get("MD_TEST_RTOL", "1e-9"))
ATOL = float(os.environ.get("MD_TEST_ATOL", "1e-14"))


def pytest_addoption(parser):
    parser.addoption("--update-golden", action="store_true",
                     help="rewrite tests/golden/*.json from the current binaries instead of comparing")


@pytest.fixture(scope="session")
def update_golden(request):
    return request.config.getoption("--update-golden")


def pytest_sessionstart(session):
    missing = [p for p in (MD2D, MD3D) if not os.access(p, os.X_OK)]
    if missing:
        pytest.exit("binaries not found (run `make`): %s" % ", ".join(missing), returncode=2)

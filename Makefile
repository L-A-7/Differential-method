# Top-level targets: build the programs, run the tests and the benchmark.
#
#   make                  build md2D and md3D
#   make test             fast test suite (~30 s): regression, physics, analytic, inputs, C units
#   make test-validation  thesis validation tables at reduced resolution (~1.5 min)
#   make update-golden    rewrite tests/golden/ from the current build (only after an intended change)
#   make bench            timing benchmark, appended to tests/bench_history.jsonl
#
# pytest is taken from the system if available, otherwise installed once into .venv/.

PYTEST_SYS := $(shell python3 -c "import pytest" 2>/dev/null && echo yes)
PYTHON := $(if $(PYTEST_SYS),python3,.venv/bin/python)
PYDEP  := $(if $(PYTEST_SYS),,.venv/bin/python)

.PHONY: all build md2D md3D test test-validation test-all update-golden bench bench-quick clean-tests

all: build

build: md2D md3D

md2D:
	$(MAKE) -C md2D

md3D:
	$(MAKE) -C md3D

.venv/bin/python:
	python3 -m venv --system-site-packages .venv
	.venv/bin/pip install --quiet pytest

test: build $(PYDEP)
	$(PYTHON) -m pytest

test-validation: build $(PYDEP)
	$(PYTHON) -m pytest -m validation -v --durations=5

test-all: build $(PYDEP)
	$(PYTHON) -m pytest -m ""

update-golden: build $(PYDEP)
	$(PYTHON) -m pytest tests/test_regression.py --update-golden -q

bench: build
	python3 tests/bench.py

bench-quick: build
	python3 tests/bench.py --quick

clean-tests:
	rm -rf .pytest_cache tests/__pycache__

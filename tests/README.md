# Test suite

Run from the repository root:

| Command | What | Time |
|---|---|---|
| `make test` | fast suite: everything below except the thesis tables | < 1 min |
| `make test-validation` | thesis validation tables at reduced resolution (`test_validation.py`) | ~1.5 min |
| `make test-all` | both | |
| `make update-golden` | rewrite `golden/` from the current build | seconds |
| `make bench` / `make bench-quick` | timing benchmark, appended to `bench_history.jsonl` | ~5 min / ~2 min |

pytest is used from the system if installed, otherwise `make` creates `.venv/` with it once.
Any pytest option works directly, e.g. `.venv/bin/python -m pytest -k md3D -x`.

## What is tested

| File | Checks |
|---|---|
| `test_regression.py` | 22 small configurations (`cases.py`) covering md2D TE/TM × RK4 / Z_INVAR / IMPROVED_RCWA × surface / stack / index map, Hamming filter, smoothing, AUTO values, near field; md3D surface, pyramid, holes, conical. Results must match `golden/*.json` to a relative 1e-9 (`MD_TEST_RTOL`). |
| `test_physics.py` | energy conservation for lossless cases, its improvement with NS, e(n) = e(−n) for a symmetric profile at normal incidence, x↔y symmetry of md3D under a 90° polarisation turn, reciprocity. |
| `test_analytic.py` | flat interface vs Fresnel (md2D and md3D, any azimuth and polarisation, dielectric and metal) and absorbing thin film vs the Airy formula, to 1e-12; RK4 convergence towards Fresnel. |
| `test_crosscode.py` | md3D with `Ny = 0` equals md2D (TE and TM, dielectric and metal); RK4 converges to Z_INVAR on a lamellar grating; Z_INVAR stays stable at large N with enough slices. |
| `test_inputs.py` | complex-number syntax, command-line overrides, clean errors (non-zero exit code) for missing or invalid keys. |
| `test_c_units.py` | C unit tests of `md_libs/md_io_utils.c` (`c/test_md_io_utils.c`). |
| `test_validation.py` | all thesis cases except the timing table (07), at reduced N or NS so that each runs in seconds; dielectrics (02, 03, 06, 09/10, 11) and metals or absorbing media (01, 04, 05, 12, 13), in 2D and 3D. Where the thesis gives a convergence table (05, 09/10, 12, 13) the comparison is at the same N; otherwise the tolerance covers the convergence gap (≈ 2–3 × the measured deviation). Full-resolution runs: `thesis_validation_cases/*/run.sh`. |

## When a test fails

- **Regression**: results changed. If the change is intended (new physics, bug fix), check that
  the physics, analytic and validation tests still pass, then `make update-golden` and commit the
  new golden files with the code change. Each golden file records the commit that produced it.
- A change of compiler flags or BLAS library may move results at the 1e-12 level; that stays far
  inside the tolerance. If it does not, investigate before relaxing `MD_TEST_RTOL`.

## Adding a test case

Add an entry to `CASES` in `cases.py` (keep it under a few seconds), run `make update-golden`,
and commit the new `golden/<name>.json`.

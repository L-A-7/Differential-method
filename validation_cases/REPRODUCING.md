# Reproducing the thesis validation cases

Each `caseXX_*/` folder holds the thesis reference values (`data.csv`, `index.html`)
and, where possible, what is needed to recompute them with the current code:

- `param.txt` - md2D/md3D parameter file (comments explain the non-obvious choices)
- `run.sh`    - regenerates any missing profile, runs the solver(s) in `caseXX/run/`,
                then prints a reproduced-vs-thesis table (`tools/compare.py`)

```sh
make -C ../md2D && make -C ../md3D        # build the solvers first
./case04_metallic_grating_vs_MFS_IM/run.sh
NMAX=15 ./case08_3d_pyramid_dielectric_multi_method/run.sh   # 3D series: choose max N
```

Requirements: the md2D/md3D build dependencies (FFTW3, BLAS, LAPACK) and python3
(standard library only).

## Status (reproduction run, September 2026)

| Case | Thesis table | Status | Agreement with thesis | Notes |
|---|---|---|---|---|
| 01 | 2.1 | reproduced (matrix column) | 8e-5 rel. (order -1) | shooting method is commented out in md2D - that column cannot be recomputed |
| 02 | 2.2 | reproduced | 1-3% on small orders, 1.5e-3 on the reciprocal pair; reciprocity itself holds to 5e-5 | echelette profile regenerated, mirrored w.r.t. old md1D `echelle.txt` |
| 03 | 2.3 | reproduced | TE: ~1e-8 abs. vs MFS at NS=5000 (~1e-6 at NS=1000); TM ~1e-6 at NS=1000 | needs the 2048-point sine (256 points alias at N=100) |
| 04 | 2.4 | reproduced | 2e-4 abs. | TM on metal needs NS >= 1000 |
| 05 | 2.5 | reproduced | RCWA column: all 6 digits; DM column ~1e-5 | RCWA column = `calcul_method = Z_INVAR`, NS=1000 |
| 06 | 3.1 | reproduced, not digit-exact | 0.01-0.3% at NS=1000 (same level as thesis DM vs GSolver) | thesis caption parameters are wrong (see `param.txt`); now computed with md3D, Ny=0 |
| 07 | 4.1 | scaling only | - | absolute times are machine-dependent |
| 08 | 4.2 | reproduced | 5e-6 to 5e-4 rel. at N=9 vs thesis N=15 | pyramid profile regenerated |
| 09 | 4.3 | reproduced N=1-9 | constant ~3.5e-4 rel. offset, same convergence | idem |
| 10 | 4.4 | reproduced N=1-9 | ~5e-4 rel. | idem |
| 11 | 4.5 | reproduced N=1-8 | 1e-8 abs. from N=6 | `md3D/pr_sin3D_256x256.txt` |
| 12 | 4.6 | reproduced N=0–2 | 1.6e-4 rel. at the same N | convergence table; higher N not run (NS=250 is costly) |
| 13 | 4.7 | reproduced N=0-12 (full table) | ~4e-5 abs. | circle index map regenerated (1024x1024) |

## Things that are not obvious from the code

- **md2D requires `fft_filter` and `delta_h`** in the parameter file, but `delta_h`
  is not used: md2D takes exactly one RK4 step per S-matrix slice, so the thesis'
  "number of integration steps" corresponds to `NS` today.
- **md2D `smoothing`** defaults to 0 (before September 2026 it was left
  uninitialised when absent, and the parser could then demand `l_smooth`).
- **Profiles are rescaled** to `[0, h]` using the `h` of the parameter file, so a
  profile file only defines the shape.
- **Normal vectors in md3D**: for `H_XY` profiles md3D computes the normal field
  itself by finite differences - no separate file is needed. For `N_XY_ZINVAR`
  index maps the normal field is currently hard-coded as radial, centred on the
  cell (`Normal_N_XY_ZINVAR` in `md3D.c` returns before the code that would read
  `norm_x`/`norm_y` from the profile file).
- **Order labels**: the thesis tables are not consistent between chapters;
  `tools/compare.py` documents the mapping to md2D/md3D order indices.
- **Speed**: 3D cost grows as N^6. With reference BLAS N=10 takes about an hour,
  N=15 about ten hours per point (pyramid, NS=20: N=5 1.5 min, N=9 22 min on a shared 8-core machine); an optimised BLAS (e.g. OpenBLAS) helps a lot.

## Where the inputs come from

The reference values are the published tables of L. Arnaud's PhD thesis
(Institut Fresnel, 2008; full text: <https://theses.hal.science/tel-00385414v1/document>).
The original parameter files and scripts of those computations are kept in the
author's archive, not in this repository. The lost profile files were
regenerated with `tools/mkprofile.py` and checked against saved outputs.

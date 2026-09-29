# Provenance and recovered history

This project's top-level git tracking had been lost. The current git history
was reconstructed in 2026 by finding and merging several older, buried
version-control traces that survived inside the working copy. This document
records what was found, what could and could not be recovered, and where the
evidence lives.

## Recovered into git history

- **Bazaar (`.bzr`) repositories** in `md2D/`, `md3D/`, `md_libs/`, `utils/`
  (and duplicated copies under the now-removed `PourAlex/`) held real,
  commit-by-commit history from 2009-2010. These were converted with
  `breezy`'s `fast-export` + `git fast-import`, relocated to their module
  subdirectory with `git-filter-repo`, and merged into this repo's history
  (see the "Merge recovered pre-GitHub history" commit). Original authors,
  dates, and messages are preserved.
- **A standalone git repository nested inside `utils/`**, remoted to
  `github.com/leoofnature/utils` (a single 2018 commit; that remote is now
  gone - "Repository not found"). Its commit was grafted onto the bzr
  `utils` history as a continuation (same author, later date). The original
  repo folder is kept on disk, untracked, at
  `utils/_old_nested_git_repo_2010-2018/`.
- **GitHub history** from `github.com/L-A-7/Differential-method` (14 commits,
  2018-2019) is preserved unchanged and merged in as well.

## Found, but NOT recoverable: CVS history

`md1D/`, `md1D/utils/`, `md1D/TESTS/`, `other/md1D_mai2007/`, and
`other/md1D_oct2006/` each contain `CVS/Root`, `CVS/Repository`, and
`CVS/Entries` files - leftover bookkeeping from an even older CVS checkout,
predating the Bazaar era.

- `CVS/Root` in all of them points to `/home/lau/cvs/` - a local
  (non-networked) CVS repository on the original machine, under a different
  username ("lau"). That path does not exist on this machine.
- No actual RCS `,v` delta files (which would hold the real revision
  content/diffs) were found anywhere in this project tree, and a broader
  search of this machine turned up none either. Only the `CVS/Entries`
  bookkeeping survived - filename, revision number, and last-sync
  timestamp per file, not the file content at each revision.
- These `CVS/*` files are tracked in git specifically because they are the
  *only* surviving trace of this era - unlike the bzr/nested-git cases
  above, there was no actual history left to import.

### What the CVS metadata tells us

Earliest dated evidence found, from `md1D/utils/CVS/Entries`:

| File | Revision | Date |
|---|---|---|
| `fft_test.c` | 1.1 | 2004-10-05 |
| `tests_old.c` | 1.1 | 2004-10-08 |
| `lire_donnees.c`, `lire_string.c` | 1.1 | 2004-10-13 |
| `lire_plot.c` | 1.1 | 2004-10-14 |
| `tests.c` | 1.1 | 2004-10-18 |
| `TESTS/param01.txt` | **1.3** | 2004-10-19 |
| `plotscript_standard` | 1.1 | 2004-10-27 |
| `script_md1D.sh` | **1.3** | 2004-10-27 |
| `plotscript_var_i` | 1.1 | 2004-10-28 |
| `profilGen.c` | **1.3** | 2004-11-15 |

Several files are already past revision 1.1 at their earliest recorded sync,
meaning the true origin predates even these October 2004 timestamps.

The core `md1D/` module (from `md1D/CVS/Entries`, and identically in
`other/md1D_mai2007/CVS/Entries` and `other/md1D_oct2006/CVS/Entries`) shows
activity from May 2005 to April 2006, with high revision numbers by the
final recorded sync (e.g. `md1D.c` at **1.25**, `Makefile` at **1.12**):

| File | Revision | Date |
|---|---|---|
| `complex.h` | 1.4 | 2005-05-24 |
| `eq_diff.c` | 1.3 | 2005-05-24 |
| `ode_solve.c` | 1.3 | 2005-07-29 |
| `md1D_in_out.h` | 1.6 | 2005-09-26 |
| `Makefile` | 1.12 | 2005-09-27 |
| `md1D_utils.c` | 1.8 | 2005-11-09 |
| `md1D_utils.h` | 1.6 | 2005-11-09 |
| `md1D.h` | 1.11 | 2005-12-07 |
| `md1D_in_out.c` | 1.21 | 2005-12-08 |
| `md1D_pilot.c` | 1.24 | 2005-12-08 |
| `md1D_pilot.h` | 1.12 | 2005-12-08 |
| `std_include.h` | 1.23 | 2005-12-08 |
| `md1D.c` | 1.25 | 2005-12-08 |
| `md1D_io_utils.c`, `md1D_io_utils.h` | 1.1 | 2006-04-24 |

### An important caveat: the three "different-vintage" folders are one snapshot

`md1D/CVS/Entries`, `other/md1D_mai2007/CVS/Entries`, and
`other/md1D_oct2006/CVS/Entries` are **byte-identical**. Despite the folder
names implying checkouts from different dates, all three were copied from
the exact same CVS sync state (nothing newer than 2006-04-24) - they were
archived/renamed by hand at later dates (Oct 2006, May 2007) without ever
re-syncing to CVS. They are not three points in the project's evolution;
they are the same point, copied three times.

### If you find an old backup

If an old backup of `/home/lau/cvs/` ever turns up (an old machine, external
drive, disk image), it would contain the actual RCS `,v` files and could be
imported the same way the `.bzr` repositories were (`cvs2git` or similar),
extending this repo's history back past 2009 to at least 2004-2005.

## RCWA/ (removed): an earlier md2D snapshot, and a lost calculation mode

The top-level `RCWA/` folder (removed after this was written) did not
implement Rigorous Coupled-Wave Analysis - despite the name, it held no
eigenvalue/modal-method code and its Makefile linked no LAPACK/BLAS. Its
`md2D.c` header identified it as the same differential method, dated
January 2007, predating the Bazaar-tracked `md2D/` history (2009-2010) and
postdating `md1D`'s CVS era (2004-2006): an earlier link in the same
lineage, using an older per-program file-naming convention (`md2D_io_utils.c`
instead of the shared `md_*` library convention) and even directly
referencing `md1D_utils.h` in `dioptre_inverse.c`.

Everything of substance in `RCWA/` was already duplicated, more completely,
in `md2D/`: the Tayeb ACES-1994 sinusoidal-aluminum-grating comparison
(`param_tayeb_sinus_metal.txt`, `res_tayeb_alu.txt`) and the SPIE-2007-linked
CD/thickness reconstruction study (`Optimization_var_lambda/`, with the same
`measured_data.txt`/`simulated_data.txt`) both exist there too, and are the
more actively maintained versions (`md2D/`'s own Tayeb comparisons were
re-verified as late as 2018). `md2D_conique.tar` was checked and found
byte-identical to the loose RCWA source - no actual conical-incidence code
despite the name.

One finding did NOT survive elsewhere: `RCWA/res_tayeb_alu.txt` and
`RCWA/sinus_alu_RCWA_*.txt` carry `calcul_method : RCWA` in their output
headers, meaning `md2D` once had an actual RCWA-algorithm mode, selectable
via a `calcul_method` parameter, used to cross-check the differential
method against a real RCWA computation for this exact benchmark. No
surviving source implements it: it's absent from `RCWA/`'s own code, and
current `md2D/md2D.c` only has a *commented-out* reference to
`calcul_method == "RCWA"` in a dead progress-bar time estimate (search
`md2D/md2D.c` for `IMPROVED_RCWA` to find it). The current
`md2D/res_tayeb_alu.txt` instead shows `calcul_method : Z_INVAR` - the
real differential method superseded and replaced that RCWA-mode
comparison at some point. So: the benchmark comparison lineage survives
in `md2D/`, but the alternate RCWA code path that once produced one side
of it is a lost capability, not just lost data.

**Correction (2026 reproduction of the thesis tables):** the capability is
not actually lost. md2D's `calcul_method = Z_INVAR` (and its alias
`IMPROVED_RCWA`) propagates each S-matrix slice as a z-invariant layer, i.e.
a staircase/RCWA computation. With `NS = 1000` it reproduces the RCWA column
of thesis Table 2.5 to all six printed digits (TE, N = 0..50; see
`thesis_validation_cases/case05_rcwa_convergence_metallic_TE/`). Only the
separate `calcul_method : RCWA` label is gone.

## mdConical/ (removed): subsumed by md3D with Ny=0

`mdConical/` (a 3-file, ~1400-line standalone program for a 1D grating at
conical incidence) was removed as redundant with `md3D/`, based on reading
both codebases (not on running them - see caveat below).

`md3D`'s `struct Param_struct` (`md3D/std_include.h`) already provisions
independent harmonic truncation in both periodic directions (`Nx`, `Ny`,
with `Lx`/`Ly`, per-order arrays `sigma_x`/`sigma_y`, and the full coupling
tensor `Qxx,Qxy,Qxz,Qyy,Qyz,Qzz`). In `md3D_pilot.c`:

```c
par->vec_size = (2*par->Nx+1)*(2*par->Ny+1);
for (i=-Ny;i<=Ny;i++)
    for (j=-Nx;j<=Nx;j++){ par->nx[...] = j; par->ny[...] = i; }
```

Setting `Ny=0` makes the outer loop run once, giving every harmonic the
same `sigma_y = sigma_y0` (independent of `Ly`) and collapsing
`vec_size` to exactly `2*Nx+1` - precisely `mdConical`'s own
`vec_size = 2*N+1`. `mdConical`'s `M_matrix` only ever uses
`Qxx, Qxz, Qyy, Qzz` (never `Qxy`/`Qyz`), and its `mdC_QMatrix` function
still carries a doc-comment copied verbatim from `md3D`'s `md3D_QMatrix`
listing `Qxy`/`Qyz` in the argument list, even though the actual function
signature had already dropped them - a leftover trace of `mdConical`
having been derived by trimming md3D's general Q-tensor routine down to
this y-invariant case. `md3D_in_out.c` reads `Ny` as a plain integer with
no lower-bound check, so `Ny=0` is not rejected.

**Caveat:** this was confirmed by reading the code, not by building both
and comparing numeric output on a shared test case. If `md3D` output with
`Ny=0` is ever needed as a stand-in for `mdConical` and the result looks
wrong, that numeric check (cheap - both had working Makefiles and test
parameter files, retrievable from this commit's parent) is the way to
verify it properly.

## md1D/ (removed): subsumed by md2D, in-plane incidence only

`md1D/` ("Methode Differentielle 1D", per its own `README.txt`) was removed
as fully subsumed by `md2D/`, confirmed by reading both codebases.

Both solve a periodic-grating diffraction problem with the same shape of
machinery (`md1D/std_include.h`'s `Param_struct` has `L`, `N`, and an
`Efficacites_struct` with per-order efficiency arrays, exactly analogous
to `md2D`'s). Two real differences were found, not just a naming change:

- **Generality, and a correction:** `md1D`'s `Param_struct` has no
  `phi_i`/`psi` fields at all - only `theta_i` and a scalar `pola` (TE or
  TM): classical, in-plane incidence, decoupled scalar TE/TM. `md2D`'s
  struct does have `phi_i`/`psi`, and `md2D_classical_FFF`'s own docstring
  claims to handle "a 1D structure with conical incidence" - but tracing
  the actual dispatch in `md2D_pilot.c` shows this is not a real, complete
  vector-coupled treatment:
  - The function-pointer selection (`if (par->pola == TE) par->M_matrix =
    M_matrix_TE; ... else par->M_matrix = M_matrix_TM;`) branches only on
    `pola`, never on `phi_i` - the grating-region calculation is always
    the scalar TE-only or TM-only matrix, never both, never coupled.
  - `psi` (which would be needed to represent a mixed/elliptical incident
    polarization for genuine conical coupling) is declared in the struct
    but never read or used anywhere in `md2D.c`/`md2D_pilot.c`.
  - `phi_i` does correctly enter the homogeneous super-/substrate
    dispersion relation via `ky_0 = sigma0*tan(phi_i)` (used in
    `kz_super`/`kz_sub`), but that same `ky_0` term is commented out in
    `md2D.c`'s grating-region field equations - so the tilt is applied
    outside the grating but ignored inside it. The `ky_0` formula itself
    carries the original author's own `/* TODO: A verifier */` (not
    verified) comment.

  So `md2D`'s conical support, as coded, is at best partial/unverified,
  not a clean superset - the actually-correct vector-coupled treatment
  for genuine conical mounting lives in `mdConical`/`md3D` (see above),
  not in `md2D`'s own `STD` path. What *is* solid: at `phi_i=0`, `ky_0`
  is identically zero either way, `psi` is irrelevant since it's unused,
  and `pola` alone selects the same scalar `M_matrix_TE`/`M_matrix_TM`
  `md1D` itself would need - so `md2D` reproduces `md1D`'s classical
  in-plane case exactly, and correctly, regardless of this caveat.
  It also means there is no calculation-cost penalty for this unused
  capability: since the dispatch never builds a larger or coupled system
  for nonzero `phi_i` in the first place, `md2D` costs exactly the same
  (same `vec_size = 2*N+1`, same scalar `M_matrix_TE`/`M_matrix_TM`)
  whether `phi_i` is zero or not - `md1D`'s removal loses no efficiency
  headroom to a "conical" code path that, in practice, isn't really there.
- **Integration algorithm:** `md1D/eq_diff.c` calls GSL's generic adaptive
  stepper (`gsl_odeiv_step_rk4`) - a general-purpose library routine.
  `md2D/md_odesolve.c` has custom analytically-derived solvers instead:
  `zinvar_P_matrix`/`zinvar_M_matrix_TM` (closed-form propagation through
  homogeneous sublayers, no numerical integration needed under
  `calcul_method = Z_INVAR`), and `implicit_rk_P_matrix` (a custom implicit
  RK exploiting the block-diagonal structure of the Psi matrices - per the
  recovered Bazaar commit history, O(N^2) instead of O(N^3), with an
  analytic rather than numerical matrix inverse). A genuine algorithmic
  upgrade, not just different tuning.

`md1D/README.txt` described basic usage: parameters in a config file
(conventionally one file per calculation type, invoked as
`md1D -param param_file_name`), with the profile defined in a separate
file. `md1D/version_a_l_endroit/Structure_md1D.txt` sketched the internal
pipeline: profile -> `k2_xyz`/`invk2_xyz` -> Fourier-transformed
(`TFk2_xyz`/`TFinvk2_xyz`) -> `S_matrix`. Both are recorded here since
they're a reasonable template for `md2D`/`md3D` usage documentation, even
though `md1D`'s own runnable files (in `TESTS/`, `tests_optim/`, and
`md1D_usage_examples/`) only made sense run against the now-removed
`md1D` binary.

Everything else in `md1D/` was confirmed redundant before removal:
`evolv/`, `version_a_l_endroit/`, `version_shakti_01/` were further
frozen duplicate snapshots of `md1D` itself (same pattern as
`PourAlex/new_DM` vs `md2D`, described above); `tests_optim/` duplicated
the CD/thickness reconstruction study already preserved, more completely,
in `md2D/Optimization_var_lambda/` (same material index files, same
vendored `liblevmar.a`); the material index files (`AIR_index.txt`,
`SI_CRISTAL_index.txt`, `POLY03_index.txt`, `OXIDE_THERM_index.txt`) were
checked and are essentially identical to the copies already in
`md2D/Optimization_var_lambda/`; and `md1D/utils/` held nothing beyond
the same generic tools already in `md2D/utils/`.

The two Octave/Matlab scripts `inverse_reflexion_TE.m` and
`inverse_reflexion_TM.m` were kept (relocated to `Applications/`) as
genuine standalone post-processing tools, not duplicated elsewhere.

`md1D/CVS/`, `md1D/utils/CVS/`, and `md1D/TESTS/CVS/` metadata is not
repeated here - it is the same CVS provenance already fully documented
above (October 2004 onward), just re-encountered at its original
location; removing it loses no information not already recorded.

### A related, unexplored lead

`~/Programs/M_files/MethodDiff/m_methodDiff/CVS` (outside this project,
in a sibling `M_files` directory) is a separate MATLAB implementation with
its own CVS tracking. It has not been investigated; it may share lineage
with this code and could be worth checking if tracing the absolute origin
further matters.

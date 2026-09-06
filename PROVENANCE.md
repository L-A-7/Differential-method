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

### A related, unexplored lead

`~/Programs/M_files/MethodDiff/m_methodDiff/CVS` (outside this project,
in a sibling `M_files` directory) is a separate MATLAB implementation with
its own CVS tracking. It has not been investigated; it may share lineage
with this code and could be worth checking if tracing the absolute origin
further matters.

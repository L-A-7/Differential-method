# Differential Method: md2D / md3D

Rigorous electromagnetic diffraction by periodic and rough structures, using the
differential method in Fourier space (S-matrix propagation, fast Fourier
factorisation with a normal-vector field, matrix Runge-Kutta integration).
The code comes from L. Arnaud's PhD work at Institut Fresnel (2008).

| Program | Structures | Illumination |
|---|---|---|
| `md2D` | invariant along y, periodic along x (surface h(x), stacks, index maps) | classical incidence, TE or TM; optional near-field maps |
| `md3D` | periodic along x and y (surface h(x,y), stacks, z-invariant index maps) | any angle, azimuth and linear polarisation; conical incidence with `Ny = 0` |

## Quick start

```sh
sudo apt install build-essential libfftw3-dev libblas-dev liblapack-dev
make -C md2D && make -C md3D
validation_cases/case04_metallic_grating_vs_MFS_IM/run.sh   # ~1 min, compares with published values
```

## Testing

```sh
make test              # ~100 fast checks, under a minute
make test-validation   # thesis validation tables at reduced resolution, ~1.5 min
make bench             # timing benchmark
```

See [`tests/README.md`](tests/README.md).

## Documentation

- [`doc/index.html`](doc/index.html): user manual ([PDF](doc/manual.pdf))
- [`doc/workflow.html`](doc/workflow.html): practical guide with recipes and troubleshooting ([PDF](doc/practical_guide.pdf))
- [`doc/optimization_report.html`](doc/optimization_report.html): what was optimized and how it was checked ([PDF](doc/optimization_report.pdf))
- [`validation_cases/review.html`](validation_cases/review.html): reproduction of the thesis validation tables with the current code
- [`validation_cases/REPRODUCING.md`](validation_cases/REPRODUCING.md): how to rerun each validation case
- [`PROVENANCE.md`](PROVENANCE.md): history of the code base
- L. Arnaud, PhD thesis, Institut Fresnel, 2008 (in French), with the full derivations and the published validation tables: <https://theses.hal.science/tel-00385414v1/document>

## Layout

| Path | Content |
|---|---|
| `md2D/`, `md3D/` | the two programs, with sample parameter and profile files |
| `md_libs/` | shared I/O, maths and utility code (symlinked into both programs) |
| `utils/` | helpers: `profilGen` (profile generator), `lire_tab` (extract arrays from results), ... |
| `validation_cases/` | 13 validation cases with inputs, run scripts and comparison tools |
| `tests/` | automated test suite (pytest) and benchmark |
| `Applications/` | a near-field example and post-processing scripts |
| `doc/`, `references/` | documentation and reference material |

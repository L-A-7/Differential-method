#!/bin/sh
# Cases 08, 09, 10: convergence series N = 1..NMAX (default 10).
# Cost grows ~N^6: N=5 ~1 min, N=10 ~1 h, N=15 ~10 h with reference BLAS.
# Usage: NMAX=15 ./run.sh
HERE=$(cd "$(dirname "$0")" && pwd); . "$HERE/../tools/common.sh"
[ -f pyramid_512x512.txt ] || python3 "$TOOLS/mkprofile.py" pyramid pyramid_512x512.txt --n 512
N=1
while [ $N -le "${NMAX:-10}" ]; do
    [ -f c08_N$N.txt ] || "$MD3D" -param param.txt -Nx $N -Ny $N -profile_name c08_N$N > c08_N$N.log 2>&1
    N=$((N+1))
done
python3 "$TOOLS/compare.py" "$HERE"
python3 "$TOOLS/compare.py" "$HERE/../case09_3d_pyramid_convergence_reflection"
python3 "$TOOLS/compare.py" "$HERE/../case10_3d_pyramid_convergence_transmission"

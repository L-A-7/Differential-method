#!/bin/sh
# Case 07: timing N = 1..NMAX (default 8), pyramid grating, NS=50.
# Absolute times depend on CPU/BLAS; compare the ~N^6 growth.
HERE=$(cd "$(dirname "$0")" && pwd); . "$HERE/../tools/common.sh"
[ -f pyramid_512x512.txt ] || python3 "$TOOLS/mkprofile.py" pyramid pyramid_512x512.txt --n 512
N=1
while [ $N -le "${NMAX:-8}" ]; do
    "$MD3D" -param param.txt -Nx $N -Ny $N -profile_name c07_N$N > c07_N$N.log 2>&1
    N=$((N+1))
done
python3 "$TOOLS/compare.py" "$HERE"

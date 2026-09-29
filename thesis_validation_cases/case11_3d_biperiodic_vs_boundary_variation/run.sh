#!/bin/sh
# Case 11: N = 1..NMAX (default 8); agreement reaches ~1e-8 by N=6.
HERE=$(cd "$(dirname "$0")" && pwd); . "$HERE/../tools/common.sh"
ln -sf "$REPO/md3D/pr_sin3D_256x256.txt" .
N=1
while [ $N -le "${NMAX:-8}" ]; do
    [ -f c11_N$N.txt ] || "$MD3D" -param param.txt -Nx $N -Ny $N -profile_name c11_N$N > c11_N$N.log 2>&1
    N=$((N+1))
done
python3 "$TOOLS/compare.py" "$HERE"

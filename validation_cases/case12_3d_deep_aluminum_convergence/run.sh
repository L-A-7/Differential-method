#!/bin/sh
# Case 12: N = 0..NMAX (default 6). NS=250 makes this ~12x slower than case 11.
HERE=$(cd "$(dirname "$0")" && pwd); . "$HERE/../tools/common.sh"
ln -sf "$REPO/md3D/pr_sin3D_256x256.txt" .
N=0
while [ $N -le "${NMAX:-6}" ]; do
    [ -f c12_N$N.txt ] || "$MD3D" -param param.txt -Nx $N -Ny $N -profile_name c12_N$N > c12_N$N.log 2>&1
    N=$((N+1))
done
python3 "$TOOLS/compare.py" "$HERE"

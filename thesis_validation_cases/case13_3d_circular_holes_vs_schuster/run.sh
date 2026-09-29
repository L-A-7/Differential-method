#!/bin/sh
# Case 13: N = 0..NMAX (default 12), ~20 min total.
HERE=$(cd "$(dirname "$0")" && pwd); . "$HERE/../tools/common.sh"
[ -f circle_R025_1024x1024.txt ] || python3 "$TOOLS/mkprofile.py" circle circle_R025_1024x1024.txt --n 1024 --R 0.25 --L 1.0
N=0
while [ $N -le "${NMAX:-12}" ]; do
    [ -f c13_N$N.txt ] || "$MD3D" -param param.txt -Nx $N -Ny $N -profile_name c13_N$N > c13_N$N.log 2>&1
    N=$((N+1))
done
python3 "$TOOLS/compare.py" "$HERE"

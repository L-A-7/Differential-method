#!/bin/sh
# Case 05: DM series (RK4, NS=5000) and RCWA series (Z_INVAR, NS=1000), ~30-60 min.
HERE=$(cd "$(dirname "$0")" && pwd); . "$HERE/../tools/common.sh"
ln -sf "$REPO/md2D/sin_2048.txt" .
NLIST="0 1 2 3 4 5 6 7 8 9 10 15 20 30 40 50"
( for N in $NLIST; do "$MD2D" -param param.txt -N $N -NS 5000 -calcul_method RK4     -profile_name c05_RK4_N$N    > /dev/null 2>&1; done ) &
( for N in $NLIST; do "$MD2D" -param param.txt -N $N -NS 1000 -calcul_method Z_INVAR -profile_name c05_ZINVAR_N$N > /dev/null 2>&1; done ) &
wait
python3 "$TOOLS/compare.py" "$HERE"

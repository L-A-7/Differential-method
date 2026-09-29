#!/bin/sh
# Case 02: ~2 min (two runs in parallel).
HERE=$(cd "$(dirname "$0")" && pwd); . "$HERE/../tools/common.sh"
python3 "$TOOLS/mkprofile.py" echelette echelette_256.txt --nx 256
"$MD2D" -param param.txt -theta_i 5.0    -profile_name c02_A > c02_A.log 2>&1 &
"$MD2D" -param param.txt -theta_i -30.60 -profile_name c02_B > c02_B.log 2>&1 &
wait
python3 "$TOOLS/compare.py" "$HERE"

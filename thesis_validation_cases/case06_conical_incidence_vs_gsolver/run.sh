#!/bin/sh
# Case 06: ~2-8 min.
HERE=$(cd "$(dirname "$0")" && pwd); . "$HERE/../tools/common.sh"
ln -sf "$REPO/md3D/carre_256_1D.txt" .
"$MD3D" -param param.txt > c06.log 2>&1
python3 "$TOOLS/compare.py" "$HERE"

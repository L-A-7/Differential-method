#!/bin/sh
# Case 01: ~1-2 min.
HERE=$(cd "$(dirname "$0")" && pwd); . "$HERE/../tools/common.sh"
ln -sf "$REPO/md2D/sin_2048.txt" .
"$MD2D" -param param.txt > c01.log 2>&1
python3 "$TOOLS/compare.py" "$HERE"

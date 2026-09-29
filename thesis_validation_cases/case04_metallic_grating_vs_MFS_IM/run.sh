#!/bin/sh
# Case 04: ~3 min (TE and TM in parallel).
HERE=$(cd "$(dirname "$0")" && pwd); . "$HERE/../tools/common.sh"
ln -sf "$REPO/md2D/sin_2048.txt" .
for P in TE TM; do "$MD2D" -param param.txt -pola $P -profile_name c04_$P > c04_$P.log 2>&1 & done
wait
python3 "$TOOLS/compare.py" "$HERE"

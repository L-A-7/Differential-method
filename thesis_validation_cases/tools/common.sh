# Sourced by each caseXX/run.sh. Sets paths and moves into caseXX/run/.
# Override the solver binaries with MD2D=... / MD3D=... if needed.
set -e
REPO=$(cd "$HERE/../.." && pwd)
TOOLS="$REPO/thesis_validation_cases/tools"
MD2D="${MD2D:-$REPO/md2D/md2D}"
MD3D="${MD3D:-$REPO/md3D/md3D}"
mkdir -p "$HERE/run"
cd "$HERE/run"
cp "$HERE/param.txt" .

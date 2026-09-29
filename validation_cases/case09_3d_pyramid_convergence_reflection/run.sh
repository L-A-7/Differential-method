#!/bin/sh
# Same runs as case 08 (results are written to case08.../run/).
HERE=$(cd "$(dirname "$0")" && pwd)
exec "$HERE/../case08_3d_pyramid_dielectric_multi_method/run.sh"

#!/bin/bash
# Installs the vendored TprParser wheel matching this build's platform and
# Python version, if one exists in conda-recipe/wheels/. TprParser provides
# optional native GROMACS .tpr support; see THIRD_PARTY_LICENSES.md for its
# license (GPLv3, separate from opengbsa's own MIT license). Its absence
# does not fail the build -- every other input format and feature works
# without it, and the caller (meta.yaml's build.script) already tolerates
# this script failing.
set -uo pipefail

PYVER=$(python -c "import sys; print(f'{sys.version_info.major}{sys.version_info.minor}')")
PLATFORM=$(python -c "import sys; print(sys.platform)")
WHEEL_DIR="$(dirname "${BASH_SOURCE[0]}")/wheels"

if [[ ! -d "$WHEEL_DIR" ]]; then
    exit 0
fi

shopt -s nocaseglob
for whl in "$WHEEL_DIR"/tprparser-*.whl; do
    [[ -e "$whl" ]] || continue

    case "$whl" in
        *cp310*) [[ "$PYVER" == "310" ]] || continue ;;
        *cp311*) [[ "$PYVER" == "311" ]] || continue ;;
        *cp312*) [[ "$PYVER" == "312" ]] || continue ;;
        *) continue ;;  # unknown Python tag
    esac

    case "$whl" in
        *manylinux*) [[ "$PLATFORM" == "linux" ]] || continue ;;
        *macosx*) [[ "$PLATFORM" == "darwin" ]] || continue ;;
        *win_amd64*) [[ "$PLATFORM" == "win32" ]] || continue ;;
        *) continue ;;  # unknown platform tag
    esac

    echo "Installing vendored TprParser wheel: $whl"
    python -m pip install --no-deps "$whl"
done
shopt -u nocaseglob

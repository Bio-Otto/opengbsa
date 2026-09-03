#!/bin/bash
# post-link script executed after opengbsa conda install
# Attempt to install vendored TprParser wheel if Python version and
# platform match. TprParser is GPLv3-licensed and vendored (not a
# declared conda/pip dependency) precisely so it stays optional: native
# .tpr support installs when a matching wheel is bundled for this
# platform/Python combination, and is silently skipped otherwise
# (opengbsa itself remains fully usable via Amber/GROMACS .top/PDB
# input either way).

PYVER=$(python -c "import sys; print(f'{sys.version_info.major}{sys.version_info.minor}')")
PLATFORM=$(python -c "import sys; print(sys.platform)")
WHEEL_DIR="$CONDA_PREFIX/share/opengbsa/tpr_wheels"
if [[ -d "$WHEEL_DIR" ]]; then
    # Glob is case-insensitive here because PyPI wheel filenames use the
    # package's declared case (TprParser-*.whl), not all-lowercase.
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

        python -m pip install --no-deps "$whl" || true
    done
    shopt -u nocaseglob
fi

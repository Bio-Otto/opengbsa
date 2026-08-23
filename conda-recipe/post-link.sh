#!/bin/bash
# post-link script executed after opengbsa conda install
# Attempt to install vendored tprparser wheel if Python version matches

PYVER=$(python -c "import sys; print(f'{sys.version_info.major}{sys.version_info.minor}')")
WHEEL_DIR="$CONDA_PREFIX/share/opengbsa/tpr_wheels"
if [[ -d "$WHEEL_DIR" ]]; then
    for whl in "$WHEEL_DIR"/tprparser-*.whl; do
        case "$whl" in
            *cp310*)
                if [[ "$PYVER" == "310" ]]; then
                    python -m pip install --no-deps "$whl" || true
                fi
                ;;
            *cp311*)
                if [[ "$PYVER" == "311" ]]; then
                    python -m pip install --no-deps "$whl" || true
                fi
                ;;
            *cp312*)
                if [[ "$PYVER" == "312" ]]; then
                    python -m pip install --no-deps "$whl" || true
                fi
                ;;
            *)
                # unknown tag
                ;;
        esac
    done
fi

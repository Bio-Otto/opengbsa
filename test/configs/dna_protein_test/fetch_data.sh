#!/usr/bin/env bash
# Downloads the IRF11-MDA5 promoter DNA complex trajectory (and the
# dataset's own gmx_MMPBSA wild-type reference output) from Zenodo record
# 14377950 ("Molecular dynamics simulation of the IRF11-DNA complex",
# CC-BY-4.0, https://doi.org/10.5281/zenodo.14377950).
set -euo pipefail
cd "$(dirname "$0")"

mkdir -p raw
BASE="https://zenodo.org/records/14377950/files"

fetch() {
    local remote_name="$1" local_name="$2"
    if [ -f "raw/${local_name}" ]; then
        echo "raw/${local_name} already exists, skipping."
        return
    fi
    echo "Downloading ${local_name}..."
    curl -L -o "raw/${local_name}" "${BASE}/${remote_name}?download=1"
}

# The archive's file names contain spaces/parens; local copies use
# underscores for easier shell handling downstream.
fetch "IRF11-MDA5promoter%20complex-traj.pdb" "IRF11-MDA5promoter_complex-traj.pdb"
fetch "gmx_MMPBSA%20output%20file%20%28wt-IRF11%29.dat" "gmx_MMPBSA_wt-IRF11.dat"

echo "Done. Raw data in raw/."
echo "Next: python3 prepare_system.py"

#!/usr/bin/env bash
# Downloads the 3dd2 protein-RNA complex topology and trajectory from
# Zenodo record 6973437 (Pavan, Bassani, Sturlese, Moro -- University of
# Padova, CC-BY-4.0), "Investigating RNA-Protein Recognition Mechanisms
# through Supervised Molecular Dynamics (SuMD) Simulations",
# https://doi.org/10.5281/zenodo.6973437.
#
# Only the 3dd2 system and its first SuMD replicate (suMD1) are fetched --
# the record also contains 3 other protein-RNA systems (4pdb, 5voe, RBD-PB6)
# and 9 more 3dd2 replicates, not needed for this example.
set -euo pipefail
cd "$(dirname "$0")"

mkdir -p raw
BASE="https://zenodo.org/records/6973437/files"

fetch() {
    local name="$1"
    if [ -f "raw/${name}" ]; then
        echo "raw/${name} already exists, skipping."
        return
    fi
    echo "Downloading ${name}..."
    curl -L -o "raw/${name}" "${BASE}/${name}?download=1"
}

fetch "3dd2.prmtop"
fetch "3dd2.pdb"
fetch "3dd2_suMD1.dcd"

echo "Done. Raw data in raw/."
echo "Next: python3 prepare_system.py"

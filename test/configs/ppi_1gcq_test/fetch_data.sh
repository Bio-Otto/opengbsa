#!/usr/bin/env bash
# Downloads the raw GROMACS complex trajectories for the 1GCQ and 2OOB
# protein-protein systems from Zenodo record 6638504 ("A dynamical view of
# protein-protein complexes: studies by molecular dynamics simulations",
# CC-BY-4.0, https://doi.org/10.5281/zenodo.6638504).
#
# Each archive is large (~1.4 GB, includes solvated + dry trajectories for
# both the bound complex and each unbound chain) since Zenodo does not offer
# partial/selective downloads of files inside a tar.gz. Only the dry
# ("without_water") complex trajectory is used by prepare_system.py --
# the rest is left in raw/ in case you also want to reproduce the
# independent Amber reference build (--build-amber-reference), which needs
# nothing extra from these archives either, but is kept intact for
# transparency/reproducibility.
set -euo pipefail
cd "$(dirname "$0")"

mkdir -p raw
DOI_BASE="https://zenodo.org/records/6638504/files"

for system in 1GCQ 2OOB; do
    archive="raw/${system}.tar.gz"
    if [ -d "raw/${system}/complex/without_water" ]; then
        echo "raw/${system}/ already extracted, skipping."
        continue
    fi
    if [ ! -f "$archive" ]; then
        echo "Downloading ${system}.tar.gz (~1.4 GB)..."
        curl -L -o "$archive" "${DOI_BASE}/${system}.tar.gz?download=1"
    fi
    echo "Extracting ${system}.tar.gz..."
    mkdir -p "raw/${system}"
    tar -xzf "$archive" -C "raw/${system}" --strip-components=1
done

echo "Done. Raw data in raw/1GCQ/ and raw/2OOB/."
echo "Next: python3 prepare_system.py --system 1gcq   (or --system 2oob)"

#!/usr/bin/env bash
# Downloads SARS-CoV-2 3CLpro + baicalein MD data from Zenodo record
# 17926575, https://doi.org/10.5281/zenodo.17926575 (CC-BY-4.0),
# "Design, Synthesis and Computational Insights of 7-Hydroxystilbene-
# Coumarin Hybrid Scaffolds as SARS-CoV-2 3CLpro Inhibitors".
#
# WARNING: the dataset's 1_MD.zip bundles all 7 ligands' MD data (7
# complexes x 300 ns trajectories) into a single 9.1 GB archive with no
# per-file selective download available from Zenodo (confirmed: the file
# server does not honor HTTP Range requests for this archive, returning
# the full object regardless). There is no way to fetch only baicalein's
# ~5 MB of topology files without downloading the entire 9.1 GB archive
# first. Budget significant time/bandwidth for this one -- easily 30-90+
# minutes depending on connection. Run it in the background
# (`nohup ./fetch_data.sh > fetch.log 2>&1 &`) if you don't want to wait
# on it interactively.
#
# Only baicalein (the smallest/representative ligand, used as this
# project's worked example) is extracted. The other 6 ligands
# (compound7a-f) follow the exact same pattern -- see README.md.
set -euo pipefail
cd "$(dirname "$0")"

mkdir -p raw
ARCHIVE="raw/1_MD.zip"

if [ -d "raw/1_MD/baicalein" ]; then
    echo "raw/1_MD/baicalein/ already extracted, skipping."
else
    if [ ! -f "$ARCHIVE" ]; then
        echo "Downloading 1_MD.zip (~9.1 GB -- this will take a while)..."
        curl -L -o "$ARCHIVE" "https://zenodo.org/records/17926575/files/1_MD.zip?download=1"
    fi
    echo "Extracting baicalein/ only..."
    unzip -o "$ARCHIVE" -d raw \
        "1_MD/baicalein/complex.prmtop" \
        "1_MD/baicalein/protein.prmtop" \
        "1_MD/baicalein/ligand.prmtop" \
        "1_MD/baicalein/dry_MD_baicalein.trj"
fi

echo "Done. Raw data in raw/1_MD/baicalein/."
echo "Next: opengbsa baicalein_config.yaml"
echo
echo "To extract additional ligands (compound7a..compound7f), re-run unzip"
echo "against raw/1_MD.zip with that ligand's subfolder name, e.g.:"
echo "  unzip -o raw/1_MD.zip -d raw '1_MD/compound7a/*'"

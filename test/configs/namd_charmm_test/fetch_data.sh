#!/usr/bin/env bash
# Downloads the CTLA4-peptide NAMD/CHARMM MD dataset from Zenodo record
# 7186684, https://doi.org/10.5281/zenodo.7186684 (CC-BY-4.0). Ships as a
# single 321 MB zip (NAMD has no split-topology/trajectory convention to
# fetch selectively, unlike the GROMACS PPI datasets' per-directory
# tarballs) -- only the topology (.psf/.pdb), CHARMM36 parameter files, and
# one trajectory (.dcd) are actually used; everything else (restart files,
# other trajectory segments, .xsc/.xst) is left in raw/ but unused.
set -euo pipefail
cd "$(dirname "$0")"

mkdir -p raw
ARCHIVE="raw/Zenodo_CTLA4_Peptide.zip"

if [ -d "raw/Zenodo_CTLA4_Peptide" ]; then
    echo "raw/Zenodo_CTLA4_Peptide/ already extracted, skipping."
else
    if [ ! -f "$ARCHIVE" ]; then
        echo "Downloading Zenodo_CTLA4_Peptide.zip (~321 MB)..."
        curl -L -o "$ARCHIVE" "https://zenodo.org/records/7186684/files/Zenodo_CTLA4_Peptide.zip?download=1"
    fi
    echo "Extracting only the files this example needs..."
    unzip -o "$ARCHIVE" -d raw \
        "Zenodo_CTLA4_Peptide/ctla4_P16_wat.psf" \
        "Zenodo_CTLA4_Peptide/ctla4_P16_wat.pdb" \
        "Zenodo_CTLA4_Peptide/par_all36m_prot.prm" \
        "Zenodo_CTLA4_Peptide/toppar_water_ions_prot.str" \
        "Zenodo_CTLA4_Peptide/output/stride50_pro_ctla4_P16.D.dcd"
fi

echo "Done. Raw data in raw/Zenodo_CTLA4_Peptide/."
echo "Next: opengbsa ctla4_peptide_config.yaml"

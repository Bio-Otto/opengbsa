#!/usr/bin/env bash
# Downloads the gp120/VRC01 antibody-antigen MD dataset from Zenodo record
# 6879091, https://doi.org/10.5281/zenodo.6879091 (CC-BY-4.0), accompanying
# Conti, Lau & Ovchinnikov, Antibodies 2022, 11(3):51
# (https://doi.org/10.3390/antib11030051).
#
# The deposit holds 30 complete systems (wild type + 29 alanine mutants of the
# VRC01 antibody), each as CHARMM .pdb + .psf + .dcd, 4,601 frames and ~1.05 GB
# per trajectory -- 31.6 GB in total. Zenodo cannot slice a .dcd server-side,
# so each trajectory must be fetched whole; we therefore fetch only the five
# systems used here, chosen to span the full experimental affinity range
# reported in Table 1 of the source paper (-6.09 to -11.78 kcal/mol):
#
#   wt       dG = -10.51   wild type
#   g54a     dG = -11.78   tightest binder in the set
#   d99a     dG = -10.45   near wild type
#   n100aa   dG =  -9.40   intermediate
#   w100ba   dG =  -6.09   weakest binder, 4.4 kcal/mol below wt
#
# System composition (measured from wt.pdb/wt.psf, not assumed):
#   19,042 atoms total; 8,827 excluding the 3,405 TIP3 waters
#   segid G120 = gp120 antigen (351 residues)
#   segid HC   = VRC01 heavy chain (111 residues)
#   segid LC   = VRC01 light chain (101 residues)
#   segid ION1/ION2 = counter-ions
# NOTE: the PDB chain-ID column is BLANK in this deposit -- receptor/ligand
# selections must use segid, not chainID.
set -euo pipefail
cd "$(dirname "$0")"

SYSTEMS=(wt g54a d99a n100aa w100ba)
BASE="https://zenodo.org/records/6879091/files"

mkdir -p raw
cd raw

for s in "${SYSTEMS[@]}"; do
    for ext in pdb psf dcd; do
        f="${s}.${ext}"
        if [ -s "$f" ]; then
            echo "$f already present, skipping."
            continue
        fi
        echo "Downloading $f ..."
        # --continue-at - resumes a partial file; the .dcd files are ~1 GB each
        curl -L --continue-at - --retry 3 --retry-delay 5 \
             -o "$f" "${BASE}/${f}?download=1"
    done
done

echo
echo "Done. Five systems in raw/ ($(du -sh . | cut -f1))."

# The .psf files need the CHARMM36m protein and water/ion parameter files, which
# the VRC01 record does not include. The CTLA4-peptide record of
# ../namd_charmm_test/ ships the same two files.
cd ..
CTLA4=../namd_charmm_test/raw/Zenodo_CTLA4_Peptide
if [ ! -s params/par_all36m_prot.prm ] || [ ! -s params/toppar_water_ions_prot.str ]; then
    if [ ! -s "$CTLA4/par_all36m_prot.prm" ]; then
        echo "Fetching the CHARMM36m parameter files with ../namd_charmm_test/fetch_data.sh ..."
        ../namd_charmm_test/fetch_data.sh
    fi
    mkdir -p params
    cp "$CTLA4/par_all36m_prot.prm" "$CTLA4/toppar_water_ions_prot.str" params/
    echo "CHARMM36m parameters copied to params/."
fi
echo "Next (from this folder): opengbsa vrc01_wt_config.yaml"

#!/usr/bin/env python3
"""
Splits the 3dd2 protein-RNA complex topology (Zenodo 6973437) into separate
receptor (protein) and ligand (RNA) prmtops for OpenGBSA's
`receptor_topology`/`ligand_topology` inputs.

Unlike the GROMACS PPI/DNA datasets, this system already ships a real
Amber .prmtop (built by the original authors via tleap, already
GB-radii-parameterized) -- no topology reconstruction is needed here, only
splitting the single complex.prmtop into its two components by atom range
(protein = atoms 0-4627, RNA = atoms 4628-5433, confirmed by this dataset's
own atom count: 4628 + 806 = 5434).

One fixup is applied after splitting: ParmEd's `Structure.strip()` zeroes
out the informational `RADIUS_SET` metadata string (leaving it as the
literal string "0") on the derived receptor/ligand prmtops, which
MMPBSA.py's own topology-consistency check rejects even though the actual
per-atom GB radii values are correctly preserved -- purely a
metadata/reporting field, not a computational one. Restored here from the
source complex.prmtop's own value.

See ../../../../OXA-MD/Publishing/rna_protein_validation/README.md for the
full validation write-up this script's output was used to produce.

Usage:
    python3 prepare_system.py
"""
import sys
from pathlib import Path

import parmed as pmd

HERE = Path(__file__).resolve().parent
RAW_DIR = HERE / "raw"
N_PROTEIN_ATOMS = 4628


def main():
    complex_prmtop = RAW_DIR / "3dd2.prmtop"
    if not complex_prmtop.exists():
        sys.exit(f"Missing {complex_prmtop} -- run ./fetch_data.sh first.")

    out_dir = HERE / "prepared"
    out_dir.mkdir(parents=True, exist_ok=True)

    complex_struct = pmd.load_file(str(complex_prmtop))
    radius_set = complex_struct.parm_data["RADIUS_SET"]
    n_total = len(complex_struct.atoms)
    print(f"Complex: {n_total} atoms "
          f"({N_PROTEIN_ATOMS} protein + {n_total - N_PROTEIN_ATOMS} RNA)")

    # Structure.strip() (in-place atom removal via an Amber-style mask) is
    # used here rather than index slicing (complex_struct[:N]), which was
    # observed to leave the saved prmtop's LENNARD_JONES_ACOEF table
    # inconsistent with its own NTYPES pointer (a ParmEd 4.3.1 slicing bug,
    # not a data problem) -- strip() correctly rebuilds all derived NB
    # tables to match the atoms actually kept.
    protein_struct = pmd.load_file(str(complex_prmtop))
    protein_struct.strip(f"@{N_PROTEIN_ATOMS + 1}-{n_total}")
    ligand_struct = pmd.load_file(str(complex_prmtop))
    ligand_struct.strip(f"@1-{N_PROTEIN_ATOMS}")

    protein_struct.parm_data["RADIUS_SET"] = radius_set
    ligand_struct.parm_data["RADIUS_SET"] = radius_set

    protein_struct.save(str(out_dir / "protein.prmtop"), overwrite=True)
    ligand_struct.save(str(out_dir / "ligand.prmtop"), overwrite=True)
    print(f"Wrote {out_dir / 'protein.prmtop'} ({len(protein_struct.atoms)} atoms)")
    print(f"Wrote {out_dir / 'ligand.prmtop'} ({len(ligand_struct.atoms)} atoms)")

    print("\nDone. Config: 3dd2_config.yaml (already points at raw/3dd2.prmtop, "
          "raw/3dd2_suMD1.dcd, and prepared/{protein,ligand}.prmtop).")


if __name__ == "__main__":
    main()

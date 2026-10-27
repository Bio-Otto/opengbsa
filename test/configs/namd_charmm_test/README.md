# PPI (NAMD/CHARMM) validation: CTLA4-peptide complex

A CTLA4-derived peptide bound to its target protein, demonstrating
OpenGBSA's **native NAMD/CHARMM Mode**: a `.psf` topology + `.pdb`
coordinates + CHARMM36 parameter files, loaded directly with no topology
reconstruction (unlike the GROMACS-sourced PPI/DNA datasets in
`../ppi_1gcq_test/` and `../dna_protein_test/`, which ship no topology at
all).

Source: Zenodo [7186684](https://doi.org/10.5281/zenodo.7186684)
(CC-BY-4.0). Ships as a single 321 MB zip containing the full NAMD run
directory (restarts, multiple trajectory segments, minimization/
equilibration output); only the production topology/coordinates and one
trajectory segment are needed here.

## Run it

```bash
cd test/configs/namd_charmm_test
./fetch_data.sh                 # downloads 321 MB, extracts ~85 MB of it
opengbsa ctla4_peptide_config.yaml
```

For a fast smoke test (3 frames instead of 300), use
`ctla4_peptide_config_test.yaml` instead -- note that with only 3 frames
from a trajectory with ~18-20 kcal/mol frame-to-frame standard deviation,
the smoke-test result will not resemble the full 300-frame mean; it only
confirms the pipeline runs end-to-end without error.

## Expected result

| | ΔH_bind (kcal/mol) |
|---|---|
| Amber MMPBSA.py (300 frames imaged with cpptraj `autoimage`) | -32.57 ± 11.41 |
| OpenGBSA (300 frames, built-in re-imaging) | -32.06 ± 11.49 |
| OpenGBSA, `amber_compatible: true` | -32.57 ± 11.41 |

Frame means ± SD. Van der Waals and electrostatic terms agree to within 0.003
kcal/mol. At default settings the polar term differs by +0.60 kcal/mol, from
the 25 Å Born-radius cutoff of Amber MMPBSA.py (`rgbmax`), and the non-polar
term by -0.08 kcal/mol, because Amber assigns the LCPO parameters of some
CHARMM atom types differently. With `amber_compatible: true` OpenGBSA applies
the same cutoff and Amber's atom typing and agrees with MMPBSA.py to within
0.004 kcal/mol.

The peptide crosses the periodic boundary in many frames of the raw
trajectory. OpenGBSA re-images the frames itself and gives the same result as
frames imaged beforehand with cpptraj.

`python test/configs/check_validation.py ctla4` (from the repository root)
runs the case and compares it with these values; `--amber-compatible`
compares with MMPBSA.py.

## What this example demonstrates

- Native Mode for CHARMM/NAMD input: `complex_pdb` pointing at a `.psf`,
  `solvated_topology` at the matching `.pdb`, and
  `forcefield_settings.charmm_params`/`charmm_coordinates` for the CHARMM36
  parameter/stream files, with no OpenMM `ForceField` parameterization.
- `binding_mode: ppi` with `chainid`-based receptor/ligand splitting for a
  peptide-protein complex (the peptide as `ligand_selection`, not a
  small-molecule `ligand_resname`).

The tutorial page of this example is
https://opengbsa.readthedocs.io/en/latest/tutorial/namd_charmm.html.

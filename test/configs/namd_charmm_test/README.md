# PPI (NAMD/CHARMM) validation: CTLA4-peptide complex

A CTLA4-derived peptide bound to its protein partner, demonstrating
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
| Amber MMPBSA.py (independent reference, 300 frames) | -17.93 ± 18.52 |
| OpenGBSA (300 frames) | -19.71 ± 20.46 |
| **Difference** | **1.78 kcal/mol** |

Full methodological write-up (including six OpenGBSA bugs found and fixed
during this validation -- all in the native CHARMM/NAMD loading path): see
this repo's sibling analysis directory,
`OXA-MD/Publishing/namd_format_validation/README.md`.

## What this example demonstrates

- Native Mode for CHARMM/NAMD input: `complex_pdb` pointing at a `.psf`,
  `solvated_topology` at the matching `.pdb`, and
  `forcefield_settings.charmm_params`/`charmm_coordinates` for the CHARMM36
  parameter/stream files -- no OpenMM `ForceField`/Coordinate-Mode
  parameterization involved.
- `binding_mode: ppi` with `chainid`-based receptor/ligand splitting for a
  peptide-protein complex (the peptide as `ligand_selection`, not a
  small-molecule `ligand_resname`).

Full methodological write-up: see this repo's sibling analysis directory,
`OXA-MD/Publishing/namd_format_validation/README.md`.

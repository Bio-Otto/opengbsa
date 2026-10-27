# PPI (protein-protein) validation: 1GCQ and 2OOB

Two independent protein-protein complexes from the same Zenodo record,
demonstrating OpenGBSA's `binding_mode: ppi` on a plain-PDB/GROMACS input
with no topology file provided (Coordinate Mode, `chainid`-based
receptor/ligand splitting instead of a small-molecule `ligand_resname`).

- **1GCQ** (Vav/GRB2 SH3 domain complex): 1984 atoms (901 + 1083)
- **2OOB** (ubiquitin/ubiquitin-ligase complex): 1895 atoms (716 + 1179)

Source: Zenodo [6638504](https://doi.org/10.5281/zenodo.6638504), "A
dynamical view of protein-protein complexes: studies by molecular dynamics
simulations" (CC-BY-4.0). Neither system ships a topology file (no
`.top`/`.itp`/`.tpr`) -- only `.gro` coordinates and an `.xtc` trajectory --
so `prepare_system.py` parameterizes both from scratch via OpenMM's
`amber14-all.xml`, same as OpenGBSA's own Coordinate Mode does internally.

## Run it

```bash
cd test/configs/ppi_1gcq_test
./fetch_data.sh                              # downloads ~2.8 GB total (both systems)
python3 prepare_system.py --system 1gcq
python3 prepare_system.py --system 2oob
opengbsa 1gcq_config.yaml
opengbsa 2oob_config.yaml
```

For a fast smoke test (3 frames instead of 300), use `1gcq_config_test.yaml`
/ `2oob_config_test.yaml` instead.

To also reproduce the independent Amber reference used to validate these
numbers (requires `tleap` on `PATH` for 1GCQ; OpenMM+ParmEd only for 2OOB).
`--gb-model` selects which of the 5 supported GB models the reference
prmtop's radii/screen set should match (default `OBC2`; see "All 5 GB
models" below):

```bash
python3 prepare_system.py --system 1gcq --build-amber-reference --gb-model OBC2
python3 prepare_system.py --system 2oob --build-amber-reference
```

## Expected result

| System | Amber MMPBSA.py ΔH_bind (kcal/mol) | OpenGBSA ΔH_bind (kcal/mol) | Difference | Difference, Amber without the 25 Å cutoff |
|---|---|---|---|---|
| 1GCQ | -37.92 ± 6.53 | -36.30 ± 6.54 | 1.62 | 0.005 |
| 2OOB | -25.47 ± 4.15 | -23.90 ± 4.17 | 1.57 | 0.005 |

300 frames each, the last 300 of a 7501-frame trajectory, OBC2 (`igb=5`);
frame means ± SD. Van der Waals, electrostatic and non-polar terms agree to
within 0.02 kcal/mol. The difference is in the polar term: Amber MMPBSA.py
computes effective Born radii only from atoms within 25 Å (`rgbmax`), OpenMM
from every atom pair. Re-running Amber's polar term without the cutoff (last
column) removes it. With `amber_compatible: true` OpenGBSA applies the same
cutoff and reproduces MMPBSA.py to within 0.002 kcal/mol.

`python test/configs/check_validation.py 1gcq 2oob` (from the repository root)
runs both systems and compares them with these values;
`--amber-compatible` compares with MMPBSA.py.

One thing to get right when reproducing the Amber side: the 1GCQ tleap
topology lists the atoms of the ACE cap in another order than the trajectory
(H1, CH3 instead of CH3, H1). `prepare_system.py --build-amber-reference`
therefore also writes `complex_prmtop_order.dcd`. Give MMPBSA.py that file,
not the original trajectory, or it evaluates the cap with two coordinates
exchanged and the electrostatic term is off by 0.065 kcal/mol.

## All five GB models (1GCQ)

Each GB model is validated against Amber MMPBSA.py on 1GCQ with a topology
that carries the radii set of the model, built with `prepare_system.py
--system 1gcq --build-amber-reference --gb-model <model>` (written to
`prepared/1gcq/amber_reference/<model>/`): `mbondi` for HCT, `mbondi2` for
OBC1 and OBC2, `bondi` for GBn and `mbondi3` for GBn2.

| GB model | Amber `igb=` | Amber ΔH_bind | OpenGBSA ΔH_bind | Difference | Difference, no cutoff |
|---|---|---|---|---|---|
| HCT | 1 | -50.33 ± 7.50 | -48.98 ± 7.51 | 1.35 | 0.004 |
| OBC1 | 2 | -38.40 ± 6.15 | -36.78 ± 6.16 | 1.62 | 0.005 |
| OBC2 | 5 | -37.92 ± 6.53 | -36.30 ± 6.54 | 1.62 | 0.005 |
| GBn | 7 | -41.73 ± 5.26 | -40.93 ± 5.26 | 0.79 | 0.003 |
| GBn2 | 8 | -42.97 ± 4.92 | -42.15 ± 4.93 | 0.82 | -0.046 |

For every model the difference is the 25 Å cutoff of the polar term in
MMPBSA.py. Without it Amber and OpenGBSA agree to within 0.005 kcal/mol, and
to within 0.05 kcal/mol for GBn2. The inputs, outputs and scripts behind these
numbers are in the data deposit of the OpenGBSA paper.

## What this example demonstrates

- `binding_mode: ppi` with `receptor_selection`/`ligand_selection` as
  `chainid N` strings, for systems where the "ligand" is a full protein
  chain rather than a small molecule.
- Coordinate Mode parameterization (no `.prmtop`/`.tpr`/`.psf` input) with
  `forcefield_settings.protein_forcefield: amber14-all.xml`, OpenMM's
  `amber14-all` force field.
- Building a usable topology from a dataset that only ships coordinates and
  a trajectory, including ambiguous histidine protonation states and
  spurious bonds from raw PDB conversion (see the docstring and comments of
  `prepare_system.py`).

The tutorial page of this example is
https://opengbsa.readthedocs.io/en/latest/tutorial/ppi_gromacs.html.

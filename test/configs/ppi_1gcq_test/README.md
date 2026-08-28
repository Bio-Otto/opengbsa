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
opengbsa run 1gcq_config.yaml
opengbsa run 2oob_config.yaml
```

For a fast smoke test (3 frames instead of 300), use `1gcq_config_test.yaml`
/ `2oob_config_test.yaml` instead.

To also reproduce the independent Amber reference used to validate these
numbers (requires `tleap` on `PATH` for 1GCQ; OpenMM+ParmEd only for 2OOB):

```bash
python3 prepare_system.py --system 1gcq --build-amber-reference
python3 prepare_system.py --system 2oob --build-amber-reference
```

## Expected result

| System | Amber reference ΔH_bind (kcal/mol) | OpenGBSA ΔH_bind (kcal/mol) | Diff |
|---|---|---|---|
| 1GCQ | -30.11 ± 6.68 | -36.12 ± 6.00 | 6.01 |
| 2OOB | -25.49 ± 4.15 | -24.46 ± 4.11 | 1.03 |

(300 frames each, last 300 of a 7501-frame trajectory.)

## What this example demonstrates

- `binding_mode: ppi` with `receptor_selection`/`ligand_selection` as
  `chainid N` strings, for systems where the "ligand" is a full protein
  chain rather than a small molecule.
- Coordinate Mode parameterization (no native `.prmtop`/`.tpr`/`.psf`
  input) via `forcefield_settings.protein_forcefield: amber14-all.xml`.
- Building a usable topology from a dataset that only ships coordinates +
  trajectory, including handling ambiguous histidine protonation states
  and avoiding spurious bond-inference artifacts from raw PDB conversion
  (see `prepare_system.py`'s docstring and inline comments for why each
  step is needed).

Full methodological write-up (including two OpenGBSA bugs found and fixed
during this validation): see this repo's sibling analysis directory,
`OXA-MD/Publishing/ppi_validation/README.md`.

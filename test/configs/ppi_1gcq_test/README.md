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
numbers (requires `tleap` on `PATH` for 1GCQ; OpenMM+ParmEd only for 2OOB).
`--gb-model` selects which of the 5 supported GB models the reference
prmtop's radii/screen set should match (default `OBC2`; see "All 5 GB
models" below):

```bash
python3 prepare_system.py --system 1gcq --build-amber-reference --gb-model OBC2
python3 prepare_system.py --system 2oob --build-amber-reference
```

## Expected result

| System | Amber reference ΔH_bind (kcal/mol) | OpenGBSA ΔH_bind (kcal/mol) | Diff |
|---|---|---|---|
| 1GCQ | -37.89 ± 6.53 | -36.12 ± 6.00 | 1.77 |
| 2OOB | -25.49 ± 4.15 | -24.46 ± 4.11 | 1.03 |

(300 frames each, last 300 of a 7501-frame trajectory. OBC2/`igb=5`.)

## All 5 GB models, validated (1GCQ)

Every GB model this project supports (`HCT`, `OBC1`, `OBC2`, `GBn`,
`GBn2`) is now independently validated against Amber MMPBSA.py on 1GCQ,
not just OBC2 -- built via `prepare_system.py --system 1gcq
--build-amber-reference --gb-model <model>` (writes to
`prepared/1gcq/amber_reference/<model>/`, each with its own correctly
radii-matched prmtop: `mbondi` for HCT, `mbondi2` for OBC1/OBC2, `bondi`
for GBn, `mbondi3` for GBn2 -- see that script's `GB_MODEL_RADII_SET`).

| GB model | Amber `igb=` | Amber ΔH_bind | OpenGBSA ΔH_bind | Diff |
|---|---|---|---|---|
| HCT | 1 | -50.31 ± 7.50 | -47.50 ± 6.72 | 2.81 |
| OBC1 | 2 | -38.38 ± 6.15 | -35.45 ± 5.57 | 2.93 |
| OBC2 | 5 | -37.89 ± 6.53 | -36.12 ± 6.00 | 1.77 |
| GBn | 7 | -41.71 ± 5.25 | -40.93 ± 5.27 | 0.78 |
| GBn2 | 8 | -42.96 ± 4.92 | -41.51 ± 4.99 | 1.45 |

(All 300 frames.) All five agree with Amber in sign and magnitude within
a tight, consistent 0.78-2.93 kcal/mol range -- comparable to or better
than the OBC2-only validation this project previously relied on.

This validation pass found and fixed four real bugs in
`mmgbsa/mmgbsa_core.py`'s GB-model handling, none of which were visible
before every model was actually exercised against an independent
reference (prior validation work only ever used `gb_model: OBC2`):

1. **`GBSACalculator._create_fallback_obc_force` only ever built an OBC2
   force**, regardless of `self.gb_model` -- HCT/OBC1/GBn/GBn2 silently
   fell back to OBC2 physics whenever this fallback path was taken
   (confirmed: all 5 `gb_model` values produced bit-identical energies on
   this exact system before the fix, since Coordinate Mode's
   `SystemGenerator` always rejects OpenMM's `implicitSolvent` kwarg and
   lands in this fallback). Fixed by dispatching to the matching
   `openmm.app.internal.customgbforces` class
   (`GBSAHCTForce`/`GBSAOBC1Force`/`GBSAOBC2Force`/`GBSAGBnForce`/
   `GBSAGBn2Force`) instead.
2. **Explicit `receptor_topology`/`ligand_topology` inputs (loaded from
   separate prmtop files) never had their GB radii synced with the
   complex's radii**, producing a large complex/receptor+ligand radii
   mismatch on any Native-Mode run using this feature (confirmed: this
   desync alone changed a real system's OBC2 result from -27.6 to -236.0
   kcal/mol on an unrelated protein-RNA validation case -- not a
   `gb_model` effect at all, just three structures disagreeing about GB
   radii on the same atoms).
3. **`GBn`'s radii set was wrong in this codebase** (`mbondi` instead of
   `bondi`) -- confirmed this isn't an OpenGBSA-specific bug: ParmEd's own
   official `createSystem(implicitSolvent=app.GBn)` path raises the exact
   same `"Radii must be between 1 and 2 Angstroms for neck lookup"` error
   when handed mbondi's (valid, Amber-standard) 0.8 A hydroxyl-hydrogen
   radii, since GBn's neck-lookup table is calibrated for bondi's radius
   range specifically.
4. **`GBn`/`GBn2` need a model-specific per-element `screen` value**
   (OpenMM's own `_SCREEN_PARAMETERS` table), not the generic
   prmtop-derived screen value every other model uses -- this was the
   single largest remaining error source, responsible for a ~13 kcal/mol
   GBn2 discrepancy against Amber before being found and fixed (GBn's
   screen values happen to be numerically closer to the generic set, so
   its error was smaller, ~3 kcal/mol, before this fix).

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

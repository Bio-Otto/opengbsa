# DNA-protein validation: IRF11-MDA5 promoter complex

A 109-residue protein (IRF11) bound to a 42-nucleotide double-stranded DNA
fragment (the MDA5 gene promoter), demonstrating OpenGBSA's
`binding_mode: ppi` on a protein-DNA complex -- the first nucleic-acid
"ligand" combination validated in this project alongside the RNA-protein
case in `../rna_protein_test/`.

Source: Zenodo [14377950](https://doi.org/10.5281/zenodo.14377950),
"Molecular dynamics simulation of the IRF11-DNA complex" (CC-BY-4.0). 3161
atoms total (1828 protein + 1333 DNA), 100 ns GROMACS MD, delivered as a
single 101-frame multi-model PDB trajectory (dry/solvent-stripped) -- no
separate topology file. The record also ships the original authors' own
`gmx_MMPBSA` reference output for the wild-type complex (fetched here too),
the only dataset in this project's validation work where a third-party
GB reference value already existed independent of anything computed here.

## Run it

```bash
cd test/configs/dna_protein_test
./fetch_data.sh                 # downloads ~24 MB
python3 prepare_system.py
opengbsa irf11_dna_config.yaml
```

For a fast smoke test (3 frames instead of 101), use
`irf11_dna_config_test.yaml` instead.

To also reproduce the independent Amber reference (OpenMM+ParmEd, no
`tleap` required):

```bash
python3 prepare_system.py --build-amber-reference
```

## Expected result

| | ΔH_bind (kcal/mol) |
|---|---|
| Dataset's own gmx_MMPBSA reference (10 frames, `raw/gmx_MMPBSA_wt-IRF11.dat`) | -127.23 ± 6.00 |
| Independent Amber reference (101 frames, this repo's own build) | -109.09 ± 26.17 |
| OpenGBSA 0.0.13 (101 frames) | -103.85 ± 26.36 |

The 5.2 kcal/mol OpenGBSA-vs-Amber difference is almost all in the polar term
(ΔEGB +6.0 of about 7200 kcal/mol): Amber MMPBSA.py computes effective Born
radii only within 25 Å (`rgbmax`, fixed by MMPBSA.py), OpenMM over every pair.
Re-running Amber's polar term without the cutoff leaves 0.58 kcal/mol (van der
Waals 0.002, electrostatics 0.25, LCPO non-polar 0.54).

Up to OpenGBSA 0.0.12 the difference was 36.7 kcal/mol, from two sources:
Coordinate Mode gave every atom one GB radius per element (every hydrogen
1.2 Å) instead of the mbondi2 set (fixed in 0.0.13; 22.5 kcal/mol here), and the
earlier Amber reference reported a non-polar term (-21.22) that a re-run of the
same recipe does not reproduce (-12.44). The Amber reference above is that
re-run; its inputs, outputs and script are in the paper's data deposit
(`coordinate_mode_references/`).

## What this example demonstrates

- `binding_mode: ppi` applied to a protein-nucleic-acid complex, using
  `chainid`-based splitting rather than a small-molecule `ligand_resname`.
- Amber's residue-template force fields (`amber14-all.xml`) natively
  support DNA residue names (`DA`/`DC`/`DG`/`DT` plus `5`/`3` terminal
  variants) with no extra configuration.
- Handling a dataset delivered as a single multi-model PDB trajectory
  (rather than a separate coordinate + trajectory file pair).
- Histidine protonation-state relabeling (HIS -> HID/HIE/HIP), needed
  because Amber-family force fields have no generic "HIS" template.

Full methodological write-up: see this repo's sibling analysis directory,
`OXA-MD/Publishing/dna_protein_validation/README.md`.

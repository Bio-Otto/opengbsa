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
| gmx_MMPBSA result shipped with the dataset (10 frames, `raw/gmx_MMPBSA_wt-IRF11.dat`) | -127.23 ± 6.00 |
| Amber MMPBSA.py (101 frames) | -109.62 ± 26.20 |
| OpenGBSA (101 frames) | -103.85 ± 26.36 |
| OpenGBSA, `amber_compatible: true` | -109.65 ± 26.20 |

OBC2, frame means ± SD. The 5.8 kcal/mol difference between OpenGBSA at
default settings and Amber is almost all in the polar term (+6.0 of about
+7,220 kcal/mol): Amber MMPBSA.py computes effective Born radii only from
atoms within 25 Å (`rgbmax`), OpenMM from every atom pair. With
`amber_compatible: true` OpenGBSA applies the same cutoff and agrees with
MMPBSA.py to -0.024 kcal/mol. That remainder comes from the Coulomb constant,
which is 3.5e-5 larger in OpenMM than in Amber: the electrostatic term
(-7,225 kcal/mol) differs by -0.25 kcal/mol, the polar term by +0.22 kcal/mol,
and the two largely cancel. The gmx_MMPBSA value of the dataset covers only 10
frames and other settings, so it is listed for reference only.

`python test/configs/check_validation.py irf11` (from the repository root)
runs the case and compares it with these values; `--amber-compatible`
compares with MMPBSA.py. The Amber reference is rebuilt with
`python3 prepare_system.py --build-amber-reference`; its inputs and outputs
are also in the data deposit of the OpenGBSA paper.

## What this example demonstrates

- `binding_mode: ppi` applied to a protein-nucleic-acid complex, using
  `chainid`-based splitting rather than a small-molecule `ligand_resname`.
- OpenMM's `amber14-all` force field natively supports DNA residue names
  (`DA`/`DC`/`DG`/`DT` plus `5`/`3` terminal variants) with no extra
  configuration.
- Handling a dataset delivered as a single multi-model PDB trajectory
  (rather than a separate coordinate + trajectory file pair).
- Histidine protonation-state relabeling (HIS -> HID/HIE/HIP), needed
  because Amber-family force fields have no generic "HIS" template.

The tutorial page of this example is
https://opengbsa.readthedocs.io/en/latest/tutorial/dna_protein.html.

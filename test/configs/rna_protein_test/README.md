# RNA-protein validation: 3dd2

A 288-residue protein bound to a 25-nucleotide RNA chain, demonstrating
OpenGBSA's handling of a protein-RNA complex via native Amber topologies
(`receptor_topology`/`ligand_topology`) rather than `binding_mode: ppi`
chain-splitting -- the complement to the DNA case in
`../dna_protein_test/`, which instead uses Coordinate Mode + `chainid`
splitting.

Source: Zenodo [6973437](https://doi.org/10.5281/zenodo.6973437) (Pavan,
Bassani, Sturlese, Moro -- University of Padova, CC-BY-4.0),
"Investigating RNA-Protein Recognition Mechanisms through Supervised
Molecular Dynamics (SuMD) Simulations". System **3dd2**: 5434 atoms total
(4628 protein + 806 RNA), a real Amber `.prmtop` already built and
GB-radii-parameterized by the original authors (no topology reconstruction
needed, unlike the GROMACS-sourced PPI/DNA datasets), plus ten independent
660-frame SuMD trajectories -- only the first (`3dd2_suMD1.dcd`) is used
here.

SuMD ("Supervised MD") trajectories simulate a binding *approach*, not an
equilibrium bound state: the protein-RNA distance narrows from ~32 A at
frame 0 to ~5-7 A by frame ~300 onward. The last 300 of 660 frames (the
close-contact/bound region) are used, matching this project's established
convention across all its validation datasets.

## Run it

```bash
cd test/configs/rna_protein_test
./fetch_data.sh                 # downloads ~44 MB (topology + 1 trajectory replicate)
python3 prepare_system.py
opengbsa 3dd2_config.yaml
```

For a fast smoke test (3 frames instead of 300), use
`3dd2_config_test.yaml` instead.

## Expected result

| | ΔH_bind (kcal/mol) |
|---|---|
| Amber MMPBSA.py (independent reference, 300 frames) | -23.49 ± 8.91 |
| OpenGBSA (300 frames) | -16.58 ± 8.97 |
| **Difference** | **6.91 kcal/mol** |

Same sign, same order of magnitude, nearly identical standard deviation --
consistent with genuine frame-to-frame conformational variance in a
still-approaching (not fully equilibrated) SuMD trajectory. Full
methodological write-up (including two OpenGBSA bugs found and fixed
during this validation, both related to nucleic-acid ligand detection):
see this repo's sibling analysis directory,
`OXA-MD/Publishing/rna_protein_validation/README.md`.

## What this example demonstrates

- Native Mode with an explicit `receptor_topology`/`ligand_topology` pair
  (rather than `binding_mode: ppi` chain splitting) for a non-small-molecule
  ligand -- here, an RNA chain.
- Amber's residue-template force fields natively support RNA residue names
  (`A`/`U`/`G`/`C` plus terminal variants) with no extra configuration.
- Splitting a single complex `.prmtop` into receptor/ligand topologies via
  ParmEd's `Structure.strip()` (an Amber-mask atom range), including a
  metadata fixup (`RADIUS_SET`) that `strip()` otherwise loses.

Full methodological write-up: see this repo's sibling analysis directory,
`OXA-MD/Publishing/rna_protein_validation/README.md`.

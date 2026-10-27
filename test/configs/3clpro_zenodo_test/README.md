# Protein-ligand validation: SARS-CoV-2 3CLpro + 7 inhibitors

Seven independent protein-ligand complexes (one protease, seven different
small-molecule inhibitors), demonstrating OpenGBSA's standard Native Mode
(Amber `.prmtop`, no topology reconstruction) on a batch of related
compounds -- baicalein plus six 7-hydroxystilbene-coumarin hybrid
scaffolds (compound7a-f).

Source: Zenodo [17926575](https://doi.org/10.5281/zenodo.17926575),
"Design, Synthesis and Computational Insights of 7-Hydroxystilbene-Coumarin
Hybrid Scaffolds as SARS-CoV-2 3CLpro Inhibitors" (CC-BY-4.0). Amber24 MD,
300 ns dry (solvent-stripped) trajectories, 30000 frames each (10 ps/frame).

**Note on the raw data size**: all 7 ligands' MD data ship bundled into a
single 9.1 GB `1_MD.zip`, with no per-file selective download available
from Zenodo for this record (confirmed: the download endpoint does not
honor HTTP `Range` requests, always returning the full object). There is
no way to fetch only one ligand's ~5 MB of topology files without
downloading the entire archive. `fetch_data.sh` extracts only `baicalein/`
(the smallest/representative case, used as this project's worked example)
but the full download is still required first -- budget real time/bandwidth
for this one (30-90+ minutes depending on connection).

## Run it

```bash
cd test/configs/3clpro_zenodo_test
./fetch_data.sh                 # downloads 9.1 GB, extracts ~3.2 GB (baicalein only)
opengbsa baicalein_config.yaml
```

To also run the other six ligands, extract their subfolders the same way
(see `fetch_data.sh`'s final message for the exact `unzip` command), then:

```bash
unzip -o raw/1_MD.zip -d raw '1_MD/compound7a/*'
opengbsa compound7a_config.yaml
# ...and so on for compound7b through compound7f
```

## Expected result

The last 300 of 30000 frames of each trajectory, OBC2, frame means ± SD in
kcal/mol. "Amber" is Amber MMPBSA.py on the same frames.

| Ligand | Amber MMPBSA.py | OpenGBSA | Difference | OpenGBSA, `amber_compatible: true` |
|---|---|---|---|---|
| baicalein | -23.89 ± 3.77 | -23.87 ± 3.79 | 0.02 | -23.888 |
| compound7a | -22.14 ± 2.46 | -22.03 ± 2.47 | 0.10 | -22.137 |
| compound7b | -30.70 ± 3.91 | -30.64 ± 3.94 | 0.06 | -30.699 |
| compound7c | -40.89 ± 4.78 | -40.92 ± 4.81 | -0.03 | -40.886 |
| compound7d | -37.88 ± 2.90 | -37.84 ± 2.91 | 0.04 | -37.879 |
| compound7e | -27.47 ± 3.04 | -27.28 ± 3.05 | 0.20 | -27.472 |
| compound7f | -46.40 ± 4.31 | -46.46 ± 4.33 | -0.07 | -46.399 |

Van der Waals, electrostatic and non-polar terms agree to within 0.001
kcal/mol. The small differences are in the polar term: Amber MMPBSA.py
computes effective Born radii only from atoms within 25 Å (`rgbmax`), OpenMM
from every atom pair. With `amber_compatible: true` OpenGBSA applies the same
cutoff and reproduces MMPBSA.py to within 0.001 kcal/mol.

`python test/configs/check_validation.py 3clpro_baicalein` (from the repository root)
runs baicalein and compares it with these values; `--amber-compatible`
compares with MMPBSA.py.
The cases `3clpro_compound7a` to `3clpro_compound7f` run the other ligands.

## What this example demonstrates

- Native Mode with an explicit `receptor_topology`/`ligand_topology` pair
  derived from a real Amber `.prmtop` (no reconstruction needed) for a
  standard small-molecule protein-ligand complex.
- Selecting a specific frame window (`frame_start`/`frame_end`) out of a
  much longer trajectory (30000 total frames, last 300 used).

One point to know when reusing this dataset with Amber: its own
`command_MMGBA.txt` passes `complex_solv.prmtop`, the solvated topology of
62451 atoms, as the `-sp` argument of MMPBSA.py. The `dry_MD_*.trj`
trajectories are already stripped of solvent and match `complex.prmtop`
(4712 atoms), so MMPBSA.py needs `complex.prmtop` there. The configurations
here and the Amber reference use it.

The tutorial page of this example is
https://opengbsa.readthedocs.io/en/latest/tutorial/protein_ligand_3clpro.html.

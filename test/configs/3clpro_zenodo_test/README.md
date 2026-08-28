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

Full-trajectory (last 300 of 30000 frames) results, from this project's
independent Amber MMPBSA.py reference vs. OpenGBSA -- see
`OXA-MD/Publishing/3CLpro_validation/README.md` for the complete 7-ligand
comparison table and per-residue decomposition analysis. Baicalein
specifically: same sign, same order of magnitude, consistent with the
other six ligands.

## What this example demonstrates

- Native Mode with an explicit `receptor_topology`/`ligand_topology` pair
  derived from a real Amber `.prmtop` (no reconstruction needed) for a
  standard small-molecule protein-ligand complex.
- Selecting a specific frame window (`frame_start`/`frame_end`) out of a
  much longer trajectory (30000 total frames, last 300 used).
- A subtlety in this specific dataset worth knowing if reusing it: its own
  `command_MMGBA.txt` passes `complex_solv.prmtop` (the solvated, 62451-atom
  topology) as Amber MMPBSA.py's `-sp` flag, which is wrong for the
  `dry_MD_*.trj` trajectories actually provided (already solvent-stripped,
  matching `complex.prmtop`'s 4712 atoms) -- this silently produces NaN
  bond/1-4 energies if not caught. This project's own configs and Amber
  reference use `complex.prmtop` correctly; not an OpenGBSA issue, just a
  documentation gap in the original dataset.

Full methodological write-up (including the one packaging-format bug found
in OpenGBSA during this validation): see this repo's sibling analysis
directory, `OXA-MD/Publishing/3CLpro_validation/README.md`.

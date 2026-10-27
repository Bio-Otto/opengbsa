# Antibody–antigen validation: HIV-1 gp120 / VRC01 Fv

An antibody bound to its antigen, exercising OpenGBSA's **native NAMD/CHARMM
Mode** (`.psf` topology + `.pdb` coordinates + CHARMM36m parameters) on the
molecular class most relevant to biologics development. Antibody–antigen
interfaces are larger and flatter than a small-molecule pocket and dominated by
polar contacts, so polar solvation carries proportionally more of the total.

Source: Zenodo [6879091](https://doi.org/10.5281/zenodo.6879091) (CC-BY-4.0),
accompanying Conti, Lau & Ovchinnikov, *Antibodies* 2022, **11**(3):51
([doi](https://doi.org/10.3390/antib11030051)). The deposit holds 30 complete
systems — wild type plus 29 alanine mutants — as `.pdb`/`.psf`/`.dcd`, 4,601
frames and ~1.05 GB per trajectory. Five are used here, chosen to span the full
experimental affinity range that paper reports.

## Run it

```bash
cd test/configs/antibody_vrc01_test
./fetch_data.sh                 # downloads 5 x 1.05 GB
python make_configs.py          # per-system configs + solvent-free Amber trajectories
python prepare_amber.py wt g54a d99a n100aa w100ba   # Amber topologies
./run_all.sh                    # both arms, all five systems, sequentially
```

`run_all.sh` is resumable: it skips any system whose output already exists.
Budget roughly 45 min per system per arm on an 8-core CPU.

## Expected result

| Variant | Amber MMPBSA.py | OpenGBSA | Difference | OpenGBSA, `amber_compatible: true` | Exp. ΔG |
|---|---|---|---|---|---|
| wt | −52.46 ± 6.78 | −36.82 ± 6.88 | 15.64 | −52.46 ± 6.78 | −10.51 |
| G54A | −60.70 ± 6.65 | −44.91 ± 6.74 | 15.79 | −60.70 ± 6.65 | −11.78 |
| D99A | −59.17 ± 7.63 | −43.60 ± 7.70 | 15.57 | −59.16 ± 7.62 | −10.45 |
| N100AA | −64.28 ± 6.95 | −48.60 ± 6.99 | 15.68 | −64.28 ± 6.95 | −9.40 |
| W100BA | −40.90 ± 6.77 | −25.43 ± 6.84 | 15.48 | −40.90 ± 6.77 | −6.09 |

All values kcal/mol, 300 frames (4301–4601), OBC2 / `igb=5`, 0.15 M salt,
frame means ± SD.

Van der Waals and electrostatic terms agree to within 0.006 kcal/mol. The
difference at default settings is almost all in the polar term (+15.7 to
+16.0 kcal/mol): Amber MMPBSA.py computes effective Born radii only from atoms
within 25 Å (`rgbmax`), OpenMM from every atom pair, and an antibody-antigen
complex has many pairs beyond 25 Å. Re-running Amber's polar term without the
cutoff leaves 0.02 kcal/mol. The non-polar term differs by −0.17 to −0.24
kcal/mol because Amber assigns the LCPO parameters of some CHARMM atom types
differently. With `amber_compatible: true` OpenGBSA applies the same cutoff
and Amber's atom typing and agrees with MMPBSA.py to within 0.005 kcal/mol.

`python test/configs/check_validation.py vrc01_wt` (from the repository root)
runs the wild type and compares it with these values; `--amber-compatible`
compares with MMPBSA.py, and `--platform CUDA` runs it on a GPU. The cases
`vrc01_g54a`, `vrc01_d99a`, `vrc01_n100aa` and `vrc01_w100ba` run the
mutants.

Experimental ΔG is from Table 1 of the source paper (RSC3/VRC01 column), shown
for context. With five variants no correlation with experiment is claimed.

## Four things about this dataset that will bite you

Each of these produces a wrong answer rather than an error, so they are worth
knowing before adapting this case to another antibody system.

1. **MDTraj cannot read the `.psf` at all.** The antibody uses Kabat numbering
   with insertion codes (`52A`, `100A`, `100B`); MDTraj parses the PSF residue
   field as a bare `int` and raises `PSFError: Could not convert residue number
   [52A]`. OpenGBSA is unaffected — it loads `.psf` through ParmEd and builds its
   MDTraj reference from the companion `.pdb`. Don't route topology loading
   through MDTraj "for consistency".

2. **The PDB chain-ID column is blank.** Chain identity lives in the segid
   column (`G120` / `HC` / `LC` / `ION1` / `ION2` / `WAT`), which MDTraj ignores
   — it reports **one** chain of 19,042 atoms. So `chainid` selections cannot
   split this system, unlike the CTLA4 case. Residue-index ranges are used
   instead, and `make_configs.py` re-derives them **per variant**: alanine
   substitution changes the antibody's atom count (3,406–3,423 across the five),
   so reusing the wild type's indices would silently mis-split the mutants.

3. **MMPBSA.py's `-sp` path fails on CHARMM residue names.** cpptraj's strip
   mask is Amber-centric (`:WAT,Cl*,CIO,…,NA,CL`) and matches none of this
   system's `TIP3`/`SOD`/`CLA`, so all 10,224 solvent atoms survive the strip and
   `sander` dies with `FATAL: NATOM mismatch in coord and topology files`.
   `make_configs.py` writes the trajectory solvent-free up front instead, and
   `-sp` is omitted entirely.

4. **A `.psf` carries no GB radii, and the CHARMM→Amber conversion invents
   none.** Every atom comes out with `solvent_radius = 0`, and `sander` with
   `igb=5` then returns **`EGB = NaN`** — while the run completes normally and
   reports every other component correctly. Only the final `DELTA TOTAL` shows
   `nan`. `prepare_amber.py` calls `changeRadii(parm, "mbondi2")` after writing
   the Amber topology; note it rejects a `CharmmPsfFile`, so write-then-reload
   is required.

## System composition

Measured from `wt.pdb`/`wt.psf`, not assumed:

| | |
|---|---|
| Atoms | 19,042 total; **8,827** excluding 3,405 TIP3 waters |
| `G120` | gp120 antigen, 351 residues, 5,398 atoms |
| `HC` | VRC01 heavy chain, 111 residues, 1,885 atoms |
| `LC` | VRC01 light chain, 101 residues, 1,535 atoms |
| Frames | 4,601 |
| Box | 12.8 nm cubic, truncated solvation shell (solute extent ~9.3 nm) |

The complex stays bound throughout: the minimum antigen–antibody heavy-atom
distance stays between 2.4 and 2.8 Å (median 2.6 Å, every 100th frame of the
4,601), so the ensemble describes a bound state rather than an association event.

CHARMM parameters are not shipped with the deposit. `params/` holds the standard
CHARMM36m protein and water/ion sets (copied from `../namd_charmm_test/`); all
45 atom types in `wt.psf` were verified to be covered by them.

## Note on parallelism

On the CPU the frames are evaluated by parallel worker processes;
`max_workers` sets their number, and these configurations set it to 4. Serial
and parallel evaluation give the same result.

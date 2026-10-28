# Normal-mode entropy reference: baicalein (3CLpro ligand)

Fixture for `test/unit/test_nmode_entropy.py`.

- `ligand.prmtop`: the baicalein ligand topology of Zenodo
  [17926575](https://doi.org/10.5281/zenodo.17926575) ("Design, Synthesis and
  Computational Insights of 7-Hydroxystilbene-Coumarin Hybrid Scaffolds as
  SARS-CoV-2 3CLpro Inhibitors", CC-BY-4.0), unchanged. GAFF, mbondi radii,
  29 atoms. It is also used by `test/configs/3clpro_zenodo_test`.
- `ligand_frame29701.pdb`: the ligand (residue UNK) of frame 29701 of
  `dry_MD_baicalein.trj`, written with cpptraj (`strip !:UNK`).
- `mmpbsa_py_nabnmode.out`: Amber's normal-mode program, as MMPBSA.py calls it
  with its default nmode settings (AmberTools 24.8):

      mmpbsa_py_nabnmode ligand_frame29701.pdb ligand.prmtop 10000 0.001 \
        'ntpr=10000, diel=C, kappa=0.000000, cut=1000, gb=1, dielc=1.000000, temp0=298.150000' \
        ligand_frame29701.nc

  Minimized energy -173.52915 kcal/mol; entropies at 298.15 K in cal/(mol K):
  translational 42.641, rotational 33.919, vibrational 49.633, total 126.194.

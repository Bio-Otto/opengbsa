# Quickstart example

OXA-family beta-lactamase 6T1H with ligand 1151, simulated for this project
with GROMACS and CHARMM36 (the topologies are in `test/data/6t1h_1151_comp/`).
`trajectory.xtc` holds 50 frames (one every 1 ns) from the last 50 ns of a
100 ns production run, with water and ions removed.

Run from the root of the repository:

    opengbsa examples/quickstart/quickstart.yaml

Expected: a mean binding enthalpy of -24.99 kcal/mol (frame-to-frame SD 4.60)
over 50 frames, the per-residue decomposition with PHE111 and ALA126 as the
strongest contributors, and an interaction-entropy estimate flagged as
unreliable (the interaction energy fluctuates by 8.8 kcal/mol). See
https://opengbsa.readthedocs.io/en/latest/quickstart.html


import mdtraj as md

pdb = 'test/3if6_test/solv_ions.pdb'
traj = md.load(pdb)
# Residues are 0-indexed in mdtraj
# The error said "residue 242". Usually OpenMM uses 0-index or 1-index depending on context.
# Let's print 241, 242, 243.

for i in range(240, 245):
    res = traj.topology.residue(i)
    print(f"Residue {i} (Index): Name={res.name}, Chain={res.chain.index}")
    for atom in res.atoms:
        print(f"  Atom: {atom.name}")

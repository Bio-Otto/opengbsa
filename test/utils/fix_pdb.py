#!/usr/bin/env python3
"""Fix PDB atom ordering to ensure residues are contiguous"""

import mdtraj as md

# Load the stripped trajectory
traj = md.load('test/complex2_stripped.pdb')

# Sort atoms by residue to ensure contiguity
topology = traj.topology

# Create a mapping of atoms sorted by residue
atom_indices = []
for residue in topology.residues:
    for atom in residue.atoms:
        atom_indices.append(atom.index)

# Reorder the trajectory
traj_reordered = traj.atom_slice(atom_indices)

# Save the fixed PDB
traj_reordered[0].save('test/complex2_fixed.pdb')
print(f"Fixed PDB saved with {traj_reordered.n_atoms} atoms")

#!/usr/bin/env python3
"""Extract clean PDB from trajectory using GRO topology"""

import mdtraj as md

# We need a topology that matches the dry trajectory (8014 atoms)
# The 6xj3.gro has 90654 atoms (solvated), but we can use it to load
# and then select only the non-water atoms

print("Loading solvated GRO...")
gro_traj = md.load('/media/bio-otto/HDD-1/Emre/MD/Complex/Dimer/ZINC000062237144/6XJ3/6xj3.gro')
print(f"GRO has {gro_traj.n_atoms} atoms")

# Select non-water, non-ion atoms
print("Selecting protein and ligand atoms...")
selection = gro_traj.topology.select('not water and not (name NA or name CL or name MG or name K or name CA)')
print(f"Selected {len(selection)} atoms")

# Slice to get just protein+ligand
stripped = gro_traj.atom_slice(selection)
print(f"Stripped topology has {stripped.n_atoms} atoms")

# Save this as reference PDB
stripped[0].save('analysis_6xj3/6xj3_clean.pdb')
print("Saved clean PDB to analysis_6xj3/6xj3_clean.pdb")

# Verify it matches the dry trajectory
print("\nVerifying against dry trajectory...")
try:
    dry_traj = md.load('analysis_6xj3/6xj3_dry.xtc', top='analysis_6xj3/6xj3_clean.pdb', frame=0)
    print(f"✓ Successfully loaded dry trajectory with clean PDB!")
    print(f"  Trajectory has {dry_traj.n_atoms} atoms, {dry_traj.n_frames} frames")
except Exception as e:
    print(f"✗ Error: {e}")

#!/usr/bin/env python3
"""
Prepare complex2 files for analysis:
1. Repair PDB topology (add hydrogens/terminals) -> complex2_fixed.pdb
2. Generate matching trajectory (subset) -> complex2_fixed.xtc
"""

import mdtraj
from openmm import app, unit
import numpy as np
import sys

# clean up imports
import warnings
warnings.filterwarnings('ignore')

input_pdb = 'test/complex2.pdb'
input_xtc = 'test/complex2_stripped.xtc'
output_pdb = 'test/complex2_fixed.pdb'
output_xtc = 'test/complex2_fixed.xtc'

# Config settings mimicking the YAML
stride = 100
max_frames = 20

print(f"Loading topology from {input_pdb}...")
try:
    pdb = app.PDBFile(input_pdb)
except Exception as e:
    print(f"Error loading PDB: {e}")
    sys.exit(1)

print(f"Loading trajectory from {input_xtc}...")
# Use MDTraj to iterate
chunk_size = 100 # Read in chunks
iter_traj = mdtraj.iterload(input_xtc, top=input_pdb, chunk=chunk_size)

# Prepare forcefield for Modeller
forcefield = app.ForceField('amber14-all.xml', 'amber14/tip3p.xml')

fixed_frames_positions = []
final_topology = None

print("Processing frames...")
frame_count = 0
total_processed = 0

for chunk in iter_traj:
    # chunk is a Trajectory object
    if frame_count >= max_frames:
        break
    
    for i in range(0, len(chunk), stride):
        if frame_count >= max_frames:
            break
            
        # Extract frame positions
        # chunk.xyz is (n_frames, n_atoms, 3) in nm
        xyz = chunk.xyz[i] * unit.nanometer
        
        # Create Modeller for this frame
        # We need to construct a fresh Modeller from the original Topology + Frame Positions
        modeller = app.Modeller(pdb.topology, xyz)
        
        # 1. Separate Ligand (UNL)
        # Assuming we want to keep UNL positions as is, but repair protein.
        # Find UNL residues
        ligand_residues = [r for r in modeller.topology.residues() if r.name == 'UNL']
        
        # Extract protein (delete ligand)
        # Note: Modeller acts in place. We should clone if we want to extract?
        # Simpler: Delete ligand, verify protein, add hydrogens, then add ligand BACK?
        # To add ligand back, we need to save the ligand atoms/positions first.
        
        # Create a separate modeller for ligand to hold it?
        # Or just extract atoms/positions.
        
        # Actually, best way:
        # Clone the modeller? `app.Modeller(modeller.topology, modeller.positions)` works.
        
        modeller_ligand = app.Modeller(modeller.topology, modeller.positions)
        protein_residues_for_del = [r for r in modeller_ligand.topology.residues() if r.name != 'UNL']
        modeller_ligand.delete(protein_residues_for_del)
        
        # Now fix protein in the main modeller
        modeller.delete(ligand_residues)
        modeller.addHydrogens(forcefield)
        
        # Now merge: Protein (Fixed) + Ligand (Original)
        # modeller now holds fixed protein
        # modeller_ligand holds original ligand
        modeller.add(modeller_ligand.topology, modeller_ligand.positions)
        
        # Capture the result
        if final_topology is None:
            final_topology = modeller.topology
            # Save the fixed PDB from the first frame
            print(f"Saving fixed topology to {output_pdb}")
            with open(output_pdb, 'w') as f:
                app.PDBFile.writeFile(final_topology, modeller.positions, f, keepIds=True)
        else:
            # Verify topology hasn't changed (atom count invariance)
            if modeller.topology.getNumAtoms() != final_topology.getNumAtoms():
                print("Warning: Topology changed between frames! This invalidates the trajectory.")
                # This could happen if addHydrogens behaves differently (e.g. histidines)
                # But usually consistent for same input.
        
        # Store positions for XTC
        # positions is list of Vec3. MDTraj needs numpy array (n_atoms, 3)
        # Convert list of Quantity to numpy
        # positions_vec = modeller.positions.value_in_unit(unit.nanometer) # OpenMM Quantity
        # Actually modeller.positions is a list of Vec3
        
        # We need to accumulate positions to write XTC later
        # Or write incrementally? MDTraj doesn't support incremental write easily without a file handle?
        # Actually it does `mdtraj.check_saver`.
        
        fixed_frames_positions.append(modeller.positions)
        frame_count += 1
        print(f"Processed frame {frame_count}")

# Convert all positions to MDTraj trajectory
if not fixed_frames_positions:
    print("No frames processed!")
    sys.exit(1)

print(f"Building final trajectory ({len(fixed_frames_positions)} frames)...")

# OpenMM positions are Quantity(List(Vec3)).
# We need (n_frames, n_atoms, 3) numpy array
n_atoms = final_topology.getNumAtoms()
xyz_array = np.zeros((len(fixed_frames_positions), n_atoms, 3))

for i, frame_pos in enumerate(fixed_frames_positions):
    # frame_pos is standard list of Vec3 (or Quantity)
    # Check type
    if isinstance(frame_pos, unit.Quantity):
        frame_pos = frame_pos.value_in_unit(unit.nanometer)
    
    # Check if numpy array
    if hasattr(frame_pos, 'shape') and hasattr(frame_pos, 'dtype'):
         xyz_array[i] = frame_pos
    else:
        for j, atom_pos in enumerate(frame_pos):
            try:
                xyz_array[i, j, :] = [atom_pos.x, atom_pos.y, atom_pos.z]
            except AttributeError:
                xyz_array[i, j, :] = atom_pos

# Create MDTraj object (we need a topology source)
# We can save PDB to disk, load it with MDTraj to get Topology object
# Or convert OpenMM topology to MDTraj topology.
# Saving PDB and loading is safest.
temp_pdb = mdtraj.load_pdb(output_pdb)
fixed_traj = mdtraj.Trajectory(xyz_array, temp_pdb.topology)

print(f"Saving trajectory to {output_xtc}...")
fixed_traj.save_xtc(output_xtc)

print("Done!")

#!/usr/bin/env python3
"""
Refine complex2_fixed.xtc by performing constrained minimization on every frame.
This relaxes the added hydrogens while keeping heavy atoms mostly fixed.
"""

from openmm import app, unit
import openmm
import mdtraj
import numpy as np

pdb_file = 'test/complex2_fixed.pdb'
xtc_file = 'test/complex2_fixed.xtc'
output_xtc = 'test/complex2_refined.xtc'

# Load Topology and Trajectory
print(f"Loading {pdb_file}...")
pdb = app.PDBFile(pdb_file)
print(f"Loading {xtc_file}...")
traj = mdtraj.load(xtc_file, top=pdb_file)

# Separate Protein and Ligand to build System
# Need to parameterize ligand again? Yes, for the system to be simulatable.
# This makes it complicated if we just want to minimize interactions.
# We can minimize Protein ONLY (strip ligand) for the trajectory if we only care about protein stability?
# But MMGBSA needs combined system. The ligand might also be clashing.

# Actually, the quickest way to robustify this:
# Use the `GBSACalculator.parameterize_protein_amber` logic which uses OpenFF for ligand.
# But that's heavy.
# Alternative: Since we know the ligand charges are gasteiger and it's small, we can just use GAFF if available?
# No, sticking to OpenFF is better.

# Let's import our `runner` code or `core` code to build the system easily?
# `core.py` has `parameterize_protein_amber` and `parameterize_ligand_openff`.
# We can reuse the `GBSACalculator` class to create the system ONCE,
# then iterate frames, set positions, minimize, save positions.

import sys
sys.path.append('.')
from mmgbsa.core import GBSACalculator

print("Initializing System builder...")
calc = GBSACalculator(charge_method='gasteiger')

# Load ligand mol
ligand_sdf = 'test/ligand2.sdf'
# Load complex pdb for structure
complex_pdb = pdb_file 

# Parameterize (create system)
print("Creating system...")
ligand_system, ligand_top, ligand_mol = calc.parameterize_ligand_openff(ligand_sdf)
protein_system, protein_top, protein_pos = calc.parameterize_protein_amber(complex_pdb, 'UNL')
combined_system = calc.create_combined_system(protein_system, ligand_system)

# Setup Integrator/Context for Minimization
# Add Restraints to Heavy Atoms
# Add Restraints
print("Adding restraints...")
restraint = openmm.CustomExternalForce("k*periodicdistance(x, y, z, x0, y0, z0)^2")
restraint.addGlobalParameter("k", 100.0*unit.kilocalories_per_mole/unit.nanometer**2) # 100 kcal/mol/nm^2
restraint.addPerParticleParameter("x0")
restraint.addPerParticleParameter("y0")
restraint.addPerParticleParameter("z0")

# Iterate over PDB topology to add restraints to heavy atoms
for atom in pdb.topology.atoms():
    # Only restrain Protein Heavy Atoms
    if atom.residue.name != 'UNL' and atom.element.name != 'hydrogen':
        pos = pdb.positions[atom.index].value_in_unit(unit.nanometer)
        restraint.addParticle(atom.index, [pos[0], pos[1], pos[2]])
    else:
        # We must add ALL particles to CustomExternalForce?
        # NO. Only constrained ones.
        pass

combined_system.addForce(restraint)

integrator = openmm.VerletIntegrator(0.001)
platform = openmm.Platform.getPlatformByName('OpenCL') # Accelerate
context = openmm.Context(combined_system, integrator, platform)

refined_positions = []

print(f"Refining {len(traj)} frames...")
for i, frame in enumerate(traj):
    xyz = frame.xyz[0] * unit.nanometer # (n_atoms, 3)
    
    # Update Restraint Reference Positions
    # We want to restrain heavy atoms to *current frame's* heavy atom positions
    # and let hydrogens relax.
    
    # Wait, CustomExternalForce: x0, y0, z0 need to be set.
    # We iterate and set params.
    
    # Only need to update params for restrained atoms (heavy).
    # This involves mapping indices.
    
    # Optimization: If we just want to fix clashes, maybe just minimize properly without restraints but specific heavy atom masses? 
    # Or just weak restraints?
    # Simple method: Set positions. Minimize.
    # If we don't restrain, the whole structure might drift or unfold if vacuum/implicit is not perfect.
    # Implicit solvent (GBSA) is active in `combined_system` (OBC2).
    # So the protein should be stable-ish.
    # But usually one restrains backbone.
    
    # New Plan: Just minimal constraints.
    # Or simplified: Just minimize with `maxIterations=100`?
    # That clears clashes but keeps structure close.
    

    # Update Restraint Reference Positions
    # Iterate through all restraints and update x0, y0, z0 from current frame positions
    # We iterate over the force's particles, knowing they correspond to heavy atoms in order.
    # But we need the index of the atom to look up position.
    # CustomExternalForce stores particle index.
    
    for k in range(restraint.getNumParticles()):
        idx, params = restraint.getParticleParameters(k) # params is list [x0, y0, z0] but initially empty/dummy
        # idx is atom index in system
        
        # Get position from xyz (numpy array, nm)
        pos = xyz[idx]
        
        # Update parameters
        restraint.setParticleParameters(k, idx, [pos[0], pos[1], pos[2]])
    
    restraint.updateParametersInContext(context)
    
    context.setPositions(xyz)
    
    openmm.LocalEnergyMinimizer.minimize(context, maxIterations=50) # Very short min
    
    state = context.getState(getPositions=True)
    refined_positions.append(state.getPositions(asNumpy=True).value_in_unit(unit.nanometer))
    
    if i % 5 == 0:
        print(f"Refined frame {i}")

# Save
print(f"Saving to {output_xtc}")
# We need topology. `pdb` variable from app.PDBFile is not good enough?
# We built `protein_top` and `ligand_top`. We need joined.
# We can use the PDB file as topology source for MDTraj since atom count/order should match.
refined_traj = mdtraj.Trajectory(np.array(refined_positions), traj.topology)
refined_traj.save_xtc(output_xtc)
print("Done.")

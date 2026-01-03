#!/usr/bin/env python3
"""
Diagnose energy breakdown for complex2.
"""

from openmm import app, unit
import openmm
import sys
import numpy as np

sys.path.append('.')
from mmgbsa.core import GBSACalculator

pdb_file = 'test/complex2_fixed.pdb'
ligand_sdf = 'test/ligand2.sdf'
ligand_resname = 'UNL'

print("Initializing Calculator...")
calc = GBSACalculator(charge_method='gasteiger')

print("Parameterizing Ligand...")
ligand_system, ligand_top, ligand_mol = calc.parameterize_ligand_openff(ligand_sdf)

print("Parameterizing Protein...")
protein_system, protein_top, protein_pos = calc.parameterize_protein_amber(pdb_file, ligand_resname)

print("Creating Combined System...")
complex_system = calc.create_combined_system(protein_system, ligand_system)

# Assign Force Groups
print("\nAssigning Force Groups...")
force_names = {}
for i, f in enumerate(complex_system.getForces()):
    f.setForceGroup(i)
    force_names[i] = f.__class__.__name__
    print(f"  Force {i}: {force_names[i]}")

# Setup Context
integrator = openmm.VerletIntegrator(0.001)
platform = openmm.Platform.getPlatformByName('OpenCL')
context = openmm.Context(complex_system, integrator, platform)

# Set Positions
print(f"Protein atoms: {protein_system.getNumParticles()}")
print(f"Ligand atoms: {ligand_system.getNumParticles()}")
print(f"Complex particles: {complex_system.getNumParticles()}")

# Get positions from PDB (assuming PDB matches joined logic)
pdb = app.PDBFile(pdb_file)
# Reconstruct combined_pos explicitly.
p_atoms = [a for a in pdb.topology.atoms() if a.residue.name != ligand_resname]
l_atoms = [a for a in pdb.topology.atoms() if a.residue.name == ligand_resname]
positions = pdb.positions
p_pos = [positions[a.index].value_in_unit(unit.nanometer) for a in p_atoms]
l_pos = [positions[a.index].value_in_unit(unit.nanometer) for a in l_atoms]
combined_pos = p_pos + l_pos

context.setPositions(combined_pos)

# Calculate Energy
state = context.getState(getEnergy=True)
total_e = state.getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole)
print(f"\nTotal Potential Energy: {total_e:.2f} kcal/mol")

print("\nEnergy Breakdown:")
total_calc = 0.0
for i in range(complex_system.getNumForces()):
    try:
        e = context.getState(getEnergy=True, groups=(1<<i)).getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole)
        print(f"Force {i} ({force_names[i]}): {e:.2f} kcal/mol")
        total_calc += e
        
        f = complex_system.getForce(i)
        if isinstance(f, openmm.NonbondedForce):
            print("  (Standard LJ + Coulomb)")
        elif isinstance(f, openmm.CustomNonbondedForce):
            print("  (Salt Screening?)")
        elif isinstance(f, openmm.CustomGBForce):
            print("  (Enhanced SA?)")
        elif isinstance(f, openmm.GBSAOBCForce):
            print("  (Polar Solvation OBC)")
            
    except Exception as ex:
        print(f"Force {i}: Error {ex}")

print(f"Summed Energy: {total_calc:.2f} kcal/mol")

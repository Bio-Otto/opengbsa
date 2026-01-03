#!/usr/bin/env python3
"""
Test creating a combined protein+ligand system for MM/GBSA
"""

import sys
sys.path.insert(0, '/home/bio-otto/Desktop/ATHENA-BACKUP/Desktop/mmpbsa/mmgbsa_v0.0.4')

from mmgbsa.core import GBSACalculator
from openff.toolkit.topology import Molecule
import openmm
from openmm import app, unit

# Create calculator
calc = GBSACalculator(
    temperature=300,
    verbose=1,
    gb_model='OBC2',
    salt_concentration=0.15,
    use_cache=True,
    protein_forcefield='charmm'
)

# Parametrize components
ligand_system, ligand_top, ligand_mol_obj = calc.parameterize_ligand_openff('analysis_6xj3/ligand.sdf')
protein_system, protein_top, protein_pos = calc.parameterize_protein_amber('analysis_6xj3/6xj3_mdtraj_clean.pdb', 'UNL')

print(f"Ligand system: {ligand_system.getNumParticles()} particles")
print(f"Protein system: {protein_system.getNumParticles()} particles")

# Create a combined system by merging forces
print("\nCreating combined system...")
combined_system = openmm.System()

# Add all particles from protein
for i in range(protein_system.getNumParticles()):
    mass = protein_system.getParticleMass(i)
    combined_system.addParticle(mass)

# Add all particles from ligand  
for i in range(ligand_system.getNumParticles()):
    mass = ligand_system.getParticleMass(i)
    combined_system.addParticle(mass)

print(f"Combined system: {combined_system.getNumParticles()} particles")

# Copy forces from protein system
for force in protein_system.getForces():
    combined_system.addForce(force)

# Copy forces from ligand system (with offset indices)
protein_particles = protein_system.getNumParticles()
for force in ligand_system.getForces():
    # This is simplified - proper implementation needs to handle force-specific index offsetting
    print(f"  Force type: {type(force).__name__}")

print(f"\n✓ Combined system created with {combined_system.getNumForces()} forces")
print("Note: This is a simplified test. Proper implementation needs force index offsetting.")

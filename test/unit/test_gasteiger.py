#!/usr/bin/env python3
"""
Test creating a combined protein+ligand system using Gasteiger charges
"""

import sys
sys.path.insert(0, '/home/bio-otto/Desktop/ATHENA-BACKUP/Desktop/mmpbsa/mmgbsa_v0.0.4')

from mmgbsa.core import GBSACalculator
from openff.toolkit.topology import Molecule
import openmm
from openmm import app, unit
import time

# Create calculator with explicit gasteiger method
print("Initializing calculator with charge_method='gasteiger'...")
calc = GBSACalculator(
    temperature=300,
    verbose=1,
    gb_model='OBC2',
    salt_concentration=0.15,
    use_cache=True,
    protein_forcefield='charmm',
    charge_method='gasteiger'
)

# Parametrize components
print("Parameterizing ligand...")
start = time.time()
ligand_system, ligand_top, ligand_mol_obj = calc.parameterize_ligand_openff('analysis_6xj3/ligand.sdf')
print(f"Ligand parameterization took {time.time()-start:.2f}s")

print("Parameterizing protein using Modeller...")
# Load complex
complex_pdb = app.PDBFile('analysis_6xj3/6xj3_mdtraj_clean.pdb')
modeller = app.Modeller(complex_pdb.topology, complex_pdb.positions)

# Delete ligand (assuming residue name UNL, or find it)
# The file probably has UNL or LIG or similar. 
# We found 'UNL' was used in previous call.
to_delete = [r for r in modeller.topology.residues() if r.name == 'UNL' or r.name == 'LIG' or r.name == 'UNK']
modeller.delete(to_delete)

# Add hydrogens
print("Adding hydrogens...")
forcefield = app.ForceField('charmm36.xml', 'charmm36/water.xml')
modeller.addHydrogens(forcefield)

# Create system
print("Creating protein system...")
protein_system = forcefield.createSystem(modeller.topology, nonbondedMethod=app.NoCutoff, constraints=app.HBonds)
protein_top = modeller.topology
protein_pos = modeller.positions

print(f"Ligand system: {ligand_system.getNumParticles()} particles")
print(f"Protein system: {protein_system.getNumParticles()} particles")

# Create a combined system by merging forces
# Note: Since we are testing, we can just use the internal create_combined_system
print("\nCreating combined system...")
complex_system = calc.create_combined_system(protein_system, ligand_system)

print(f"\n✓ Complex system created successfully with {complex_system.getNumParticles()} particles")

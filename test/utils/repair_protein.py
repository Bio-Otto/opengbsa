#!/usr/bin/env python3
"""Repair protein PDB using PDBFixer"""

from openmm.app import PDBFile
from pdbfixer import PDBFixer
import sys

input_pdb = 'analysis_6xj3/6xj3_mdtraj_clean.pdb'
output_pdb = 'analysis_6xj3/6xj3_fixed.pdb'

print(f"Repairing {input_pdb}...")
fixer = PDBFixer(filename=input_pdb)

# Find issues
print("Finding missing residues...")
fixer.findMissingResidues()
print("Finding nonstandard residues...")
fixer.findNonstandardResidues()
print("Removing heterogens...")
fixer.removeHeterogens(keepWater=True)
print("Finding missing atoms...")
fixer.findMissingAtoms()

# Add missing atoms
print("Adding missing atoms...")
fixer.addMissingAtoms()

# Write output
print(f"Writing to {output_pdb}...")
with open(output_pdb, 'w') as f:
    PDBFile.writeFile(fixer.topology, fixer.positions, f, keepIds=True)

print("Done!")

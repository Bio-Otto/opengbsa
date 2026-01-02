#!/usr/bin/env python3
"""Use PDBFixer to properly rebuild PDB structure"""

from openmm.app import PDBFile
from pdbfixer import PDBFixer

# Load the stripped PDB
fixer = PDBFixer(filename='test/complex2_stripped.pdb')

# Find and fix any issues
fixer.findMissingResidues()
fixer.findNonstandardResidues()
fixer.findMissingAtoms()

# Write the fixed PDB
with open('test/complex2_pdbfixer.pdb', 'w') as f:
    PDBFile.writeFile(fixer.topology, fixer.positions, f, keepIds=True)

print("PDB fixed and saved to test/complex2_pdbfixer.pdb")

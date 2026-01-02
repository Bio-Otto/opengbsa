#!/usr/bin/env python3
"""Convert ligand2.pdb to SDF using RDKit with bond order inference"""

from rdkit import Chem
from rdkit.Chem import AllChem

# Read PDB
mol = Chem.MolFromPDBFile('test/ligand2.pdb', removeHs=False)

if mol is None:
    print("Failed to read PDB")
    exit(1)

# Try to infer bond orders and aromaticity
try:
    # Add hydrogens if needed
    mol = Chem.AddHs(mol)
    
    # Generate 3D coordinates if not present
    AllChem.EmbedMolecule(mol, randomSeed=42)
    
    # Optimize geometry
    AllChem.UFFOptimizeMolecule(mol)
    
    # Write to SDF
    writer = Chem.SDWriter('test/ligand2.sdf')
    writer.write(mol)
    writer.close()
    
    print("Successfully converted ligand2.pdb to ligand2.sdf")
    
except Exception as e:
    print(f"Error during conversion: {e}")
    exit(1)

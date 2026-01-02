from rdkit import Chem
import sys
import os

pdb_file = "ligand_ac.pdb"
sdf_file = "ligand.sdf"

if not os.path.exists(pdb_file):
    print(f"Error: {pdb_file} not found")
if not os.path.exists(sdf_file):
    pass # To avoid issues if file exists

mol = Chem.MolFromPDBFile(pdb_file, removeHs=False, sanitize=False)
if mol is None:
    mol = Chem.MolFromPDBFile(pdb_file, removeHs=False)
if mol is None:
    print("Error: Could not read PDB file with RDKit")
    sys.exit(1)

writer = Chem.SDWriter(sdf_file)
writer.SetKekulize(False) # Attempt to write even if aromaticity is complex
try:
    writer.write(mol)
    print(f"Successfully converted {pdb_file} to {sdf_file}")
except Exception as e:
    print(f"Error writing SDF: {e}")
finally:
    writer.close()

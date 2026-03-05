
from rdkit import Chem
from rdkit.Chem import AllChem
import sys
import os

input_pdb = 'test/3if6_test/ligand.pdb'
output_sdf = 'test/3if6_test/ligand.sdf'

print(f"Converting {input_pdb} to {output_sdf}...")

mol = Chem.MolFromPDBFile(input_pdb, removeHs=False)
if mol is None:
    print("Error: Could not read PDB file.")
    sys.exit(1)

try:
    # Attempt to assign bond orders from template if available, or just use 3D conformer
    # Assign stereochemistry from 3D coordinates
    Chem.AssignStereochemistryFrom3D(mol)
    Chem.AssignAtomChiralTagsFromStructure(mol)
    
    # Check centers
    centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    print(f"Chiral centers: {centers}")
    
    # Write to SDF with V3000 to ensure stereochemistry is preserved? No, V2000 is safer.
    writer = Chem.SDWriter(output_sdf)
    writer.SetKekulize(False) # Try to keep original bond orders?
    writer.write(mol)
    writer.close()
    print("Success: SDF created.")
except Exception as e:
    print(f"Conversion failed: {e}")
    sys.exit(1)

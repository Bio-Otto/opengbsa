import parmed as pmd
import sys
import os

input_file = "ligand.mol2"
output_file = "ligand.sdf"

if not os.path.exists(input_file):
    print(f"Error: {input_file} not found")
    sys.exit(1)

try:
    mol = pmd.load_file(input_file)
    mol.save(output_file)
    print(f"Successfully converted {input_file} to {output_file}")
except Exception as e:
    print(f"Error with ParmEd: {e}")
    sys.exit(1)

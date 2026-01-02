from openff.toolkit.topology import Molecule
import sys

ligand_mol = "ligand.sdf"
try:
    mol = Molecule.from_file(ligand_mol, file_format='sdf', allow_undefined_stereo=True)
    # OpenFF 'from_file' returns a list if multiple are found, or a Molecule?
    # Actually current toolkit returns a list of Molecules if using RDKit backend usually.
    # But usually core.py handles it? Wait, core.py assumes it returns a Molecule object directly?
    # Let's check type.
    if isinstance(mol, list):
        print(f"Loaded list of {len(mol)} molecules")
        mol = mol[0]
    
    print(f"Atom count: {mol.n_atoms}")
except Exception as e:
    print(f"Error: {e}")


from openff.toolkit.topology import Molecule
import numpy as np

def check_ligand_charge(sdf_file):
    print(f"Loading ligand from {sdf_file}")
    mol = Molecule.from_file(sdf_file)
    
    # It might load multiple conformers/mols, take first
    if isinstance(mol, list):
        mol = mol[0]
        
    print(f"Name: {mol.name}")
    print(f"Atoms: {mol.n_atoms}")
    
    charges = mol.partial_charges
    if charges is None:
        print("No partial charges found in SDF!")
        return

    # charges is a Quantity array
    q_vals = charges.magnitude # assuming elementary charge units likely
    # Usually OpenFF uses elementary charge e
    
    net_charge = np.sum(q_vals)
    print(f"Net Charge: {net_charge:.4f} e")
    
    min_q = np.min(q_vals)
    max_q = np.max(q_vals)
    print(f"Min Partial Charge: {min_q:.4f}")
    print(f"Max Partial Charge: {max_q:.4f}")
    
    # Check for extreme charges
    if min_q < -1.0 or max_q > 1.0:
        print("Warning: Some atoms have high partial charge magnitude (>1.0).")
        
    # Print high charge atoms
    for i, q in enumerate(q_vals):
        if abs(q) > 0.8:
            print(f"  Atom {i} ({mol.atoms[i].symbol}): {q:.4f}")

if __name__ == "__main__":
    import sys
    check_ligand_charge(sys.argv[1])


import mdtraj as md
import sys

input_gro = 'test/3if6_test/solv_ions.gro'
output_pdb = 'test/3if6_test/solv_ions.pdb'

print(f"Loading {input_gro}...")
try:
    traj = md.load(input_gro)
    print(f"Loaded. Atoms: {traj.n_atoms}, Residues: {traj.n_residues}")
    
    # Save as PDB
    traj.save(output_pdb)
    print(f"Converted to {output_pdb}")
    
except Exception as e:
    print(f"Error: {e}")
    sys.exit(1)

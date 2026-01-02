
import mdtraj as md
import sys

input_pdb = 'test/3if6_test/3if6.pdb'
output_pdb = 'test/3if6_test/3if6_fixed.pdb'

print(f"Loading {input_pdb} with MDTraj...")
try:
    traj = md.load(input_pdb)
    print(f"Loaded. Atoms: {traj.n_atoms}, Residues: {traj.n_residues}")
    
    # Save to new PDB
    # MDTraj writes standard PDB format, often fixing connectivity issues implicitly
    # by not writing weird CONECT records, or writing standard ones.
    traj.save(output_pdb)
    print(f"Saved to {output_pdb}")
    
    # Validate
    traj2 = md.load(output_pdb)
    print(f"Reloaded. Atoms: {traj2.n_atoms}")
    
    if traj.n_atoms == traj2.n_atoms:
        print("Success: Atom count preserved.")
        # Check standard OpenMM loading (the ultimate test)
        from openmm import app
        try:
            pdb = app.PDBFile(output_pdb)
            print("Success: OpenMM PDBFile can read the fixed file!")
        except Exception as e:
            print(f"Failure: OpenMM still rejects it: {e}")
            sys.exit(1)
    else:
        print("Failure: Atom count changed!")
        sys.exit(1)

except Exception as e:
    print(f"Error: {e}")
    sys.exit(1)

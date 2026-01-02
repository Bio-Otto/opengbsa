import mdtraj as md
import os

def strip_trajectory():
    print("Loading trajectory with PDB topology...")
    # Load solvated trajectory with PDB topology
    try:
        t = md.load('test/complex2.xtc', top='test/complex2.pdb')
        print(f"Loaded trajectory: {t.n_frames} frames, {t.n_atoms} atoms")
    except Exception as e:
        print(f"Error loading trajectory: {e}")
        return

    print("Identifying atoms to keep...")
    # Select protein and ligand
    # Exclude water (HOH, SOL, TIP3) and ions (NA, CL, etc)
    selection = t.topology.select('not water and not (name NA or name CL or name MG or name K)')
    
    print(f"Selected {len(selection)} atoms to keep")
    
    print("Slicing trajectory...")
    stripped_t = t.atom_slice(selection)
    
    output_traj = 'test/complex2_stripped.xtc'
    output_pdb = 'test/complex2_stripped.pdb'
    
    print(f"Saving stripped trajectory to {output_traj}...")
    stripped_t.save(output_traj)
    
    print(f"Saving stripped topology to {output_pdb}...")
    stripped_t[0].save(output_pdb)
    
    print("Done!")

if __name__ == "__main__":
    strip_trajectory()

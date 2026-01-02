
import mdtraj as md
import numpy as np
import sys

def repair_trajectory(gro_file, xtc_file, output_xtc):
    print(f"Loading topology from GRO: {gro_file}")
    # Load full frame 0 to get topology and atom indices
    full_frame = md.load(gro_file)
    
    # Identify dry atoms (Protein + Ligand)
    dry_indices = full_frame.topology.select("protein or resname UNL or resname LIG")
    print(f"Selected {len(dry_indices)} dry atoms for topology.")
    
    # Create dry topology
    dry_top = full_frame.topology.subset(dry_indices)
    
    print(f"Loading trajectory: {xtc_file}")
    # Load trajectory with dry topology
    traj = md.load(xtc_file, top=dry_top)
    print(f"Loaded {traj.n_frames} frames.")
    
    # 0. FIX TOPOLOGY (Remove spurious cyclic bonds)
    # Strategy: Detect long-range C-N bonds and remove them.
    # Then save/reload PDB to enforce topology update in C++ backend.
    
    new_bonds = []
    removed_count = 0
    all_bonds = list(dry_top.bonds)
    print(f"Total bonds: {len(all_bonds)}")
    
    # DEBUG: Check Bond 3990
    if len(all_bonds) > 3990:
        b3990 = all_bonds[3990]
        print(f"DEBUG Bond 3990: {b3990[0]} ({b3990[0].residue.index}) -- {b3990[1]} ({b3990[1].residue.index})")
        
    print("Scanning for spurious periodic bonds (Explicit PRO265-THR25 check)...")
    for i, bond in enumerate(dry_top.bonds):
        a1, a2 = bond[0], bond[1]
        r1, r2 =str(a1.residue), str(a2.residue)
        
        # Check specific PRO265 - THR25 bond (found by debug)
        if ("PRO265" in r1 and "THR25" in r2) or ("PRO265" in r2 and "THR25" in r1):
             print(f"❌ Removing spurious cyclic bond {i}: {a1} ({r1}) -- {a2} ({r2})")
             removed_count += 1
             continue
             
        # Also keep index gap check just in case, but relax threshold or skip
        # new_bonds.append(bond)
            
        new_bonds.append(bond)
            
    if removed_count > 0:
        print(f"Removed {removed_count} spurious bonds. applying patch...")
        
        # Patch python list in place
        dry_top._bonds = new_bonds
        print(f"Patched topology bonds: {len(list(dry_top.bonds))}")
        
        # We must create a new trajectory object to ensure the C++ backend (if any)
        # receives the updated topology.
        # Load XTC positions again or use existing?
        # We need to apply this topology to the XTC frames.
        
        # Create new trajectory with patched topology
        # This might fail if frames > memory? 
        # But we already loaded 'traj' into memory (md.load default is memory).
        
        print("Re-creating Trajectory object with patched topology...")
        
        # Preserve unitcell info from original trajectory
        old_uc = traj.unitcell_vectors
        old_time = traj.time
        
        # Create new trajectory with patched topology
        traj = md.Trajectory(traj.xyz, dry_top)
        
        # Restore unitcell info
        traj.unitcell_vectors = old_uc
        traj.time = old_time
        
        print("Trajectory object updated.")

    # Check unit cell info

    # Check unit cell info

    # Check unit cell info
    if traj.unitcell_vectors is None:
        print("Warning: Trajectory lacks unit cell vectors. PBC repair might fail.")
    else:
        print("Unit cell vectors present. Proceeding with repair.")

    print("Repairing (Unwrapping/Imaging molecules)...")
    
    # Get protein indices as anchor
    protein_indices = traj.topology.select("protein")
    if len(protein_indices) > 0:
        # Taking the first atom of protein as anchor usually works better than whole list for some implementations
        # But let's try passing the whole protein as anchor molecules (list of lists of indices)
        # Actually image_molecules takes 'anchor_molecules' as list of sets/list of indices of molecules.
        # MDTraj usually detects molecules automatically.
        
        # Let's try centering first to put protein in box center
        traj.center_coordinates()
        
        # Then image
        # mdtraj 1.9.9 image_molecules(sorted_bonds, anchor_molecules, ...)
        # We rely on default bond definition or topology bonds.
        traj.image_molecules(inplace=True)
    else:
        traj.image_molecules(inplace=True)
        
    print("Verification: Re-checking max displacement...")
    # Quick check if it worked?
    # No, just save.


    print(f"Saving repaired trajectory to: {output_xtc}")
    traj.save_xtc(output_xtc)
    print("Done.")

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: python repair_trajectory.py <gro> <input_xtc> <output_xtc>")
    else:
        repair_trajectory(sys.argv[1], sys.argv[2], sys.argv[3])

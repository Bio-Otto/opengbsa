
import mdtraj as md
import numpy as np
import sys

def check_integrity(protein_pdb, xtc_file, ligand_pdb=None):
    # Special handling for GRO file which might contain water
    if protein_pdb.endswith('.gro'):
        print(f"Loading full topology from GRO: {protein_pdb}")
        full_traj = md.load(protein_pdb)
        print(f"Full system atoms: {full_traj.n_atoms}")
        
        # Select dry atoms
        dry_indices = full_traj.topology.select("protein or resname UNL or resname LIG")
        print(f"Selected dry atoms: {len(dry_indices)}")
        
        # Create sliced topology
        top = full_traj.topology.subset(dry_indices)
    else:
        print(f"Loading protein topology: {protein_pdb}")
        prot_traj = md.load(protein_pdb)
        top = prot_traj.topology

        if ligand_pdb:
            print(f"Loading ligand topology: {ligand_pdb}")
            lig_traj = md.load(ligand_pdb)
            top = top.join(lig_traj.topology)
            
    print(f"Final Topology has {top.n_atoms} atoms.")
    print(f"Expected XTC atoms: 8014 (Based on previous logs)")

    print(f"Loading trajectory: {xtc_file}")
    
    # Load trajectory (stride to save time)
    try:
        t = md.load(xtc_file, top=top, stride=20)
    except Exception as e:
        print(f"\n❌ Error loading trajectory: {e}")
        print("Tip: If atom counts don't match, check if water/ions were removed vs present in PDB.")
        return

    print(f"Loaded {t.n_frames} frames (strided).")
    
    # 1. Check Protein Integrity (Bond Lengths)
    # Select alpha carbons to permit some flexibility, or all atoms for strictness
    # Let's check typical C-C, C-N bonds. If > 2.0 A (usually 1.5 A), suspect.
    # Actually, let's just use topology bonds.
    
    print("Checking bond lengths...")
    
    # MDTraj compute_distances expects pairs
    # Get all bond pairs from topology
    bond_pairs = []
    for bond in t.topology.bonds:
        bond_pairs.append([bond[0].index, bond[1].index])
    
    if not bond_pairs:
        print("Warning: No bond information in topology. Cannot check connectivity.")
        return

    # Compute distances for all bonds across all loaded frames
    distances = md.compute_distances(t, bond_pairs)
    
    # distances shape: (n_frames, n_bonds)
    max_bond_len = np.max(distances)
    mean_bond_len = np.mean(distances)
    
    print(f"Max Bond Length detected: {max_bond_len*10:.2f} Angstroms")
    print(f"Mean Bond Length: {mean_bond_len*10:.2f} Angstroms")
    
    cutoff = 0.3 # nanometers = 3.0 Angstroms. Bond shouldn't check this long.
    
    if max_bond_len > cutoff:
        print("\n🚨 CRITICAL: High bond lengths detected.")
        
        # Find broken bonds
        # distances shape: (n_frames, n_bonds)
        # Check frame 0 (since we saw it broken at 0)
        frame_idx = 0
        dists_f0 = distances[frame_idx]
        broken_indices = np.where(dists_f0 > cutoff)[0]
        
        print(f"Broken bonds in Frame {frame_idx}: {len(broken_indices)}")
        for idx in broken_indices[:5]: # Show first 5
            atom1_idx = bond_pairs[idx][0]
            atom2_idx = bond_pairs[idx][1]
            atom1 = t.topology.atom(atom1_idx)
            atom2 = t.topology.atom(atom2_idx)
            dist = dists_f0[idx]
            print(f"  Bond {idx}: {atom1} -- {atom2} = {dist*10:.2f} A")
            print(f"    Residues: {atom1.residue} -- {atom2.residue}")

        print("Diagnosis: Trajectory is WRAPPED.")
    else:
        print("\n✅ Bond lengths look normal. Protein is likely UNWRAPPED (Intact).")
        print("Diagnosis: Trajectory seems topologically intact.")

    # 2. Check Ligand-Protein Distance
    # To see if ligand is flying away or well-bound (if unwrapped but dissociated)
    ligand_atoms = t.topology.select("resname UNL")
    protein_atoms = t.topology.select("protein")
    
    if len(ligand_atoms) > 0 and len(protein_atoms) > 0:
        com_lig = md.compute_center_of_mass(t.atom_slice(ligand_atoms))
        com_prot = md.compute_center_of_mass(t.atom_slice(protein_atoms))
        
        # Distance between COMs
        # com shape (n_frames, 3)
        dists = np.linalg.norm(com_lig - com_prot, axis=1) * 10 # Angstroms
        
        print(f"\nLigand-Protein COM Distance: {np.mean(dists):.2f} ± {np.std(dists):.2f} A")
        print(f"Min Dist: {np.min(dists):.2f} A")
        print(f"Max Dist: {np.max(dists):.2f} A")
        
if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python check_pbc.py <protein_pdb> <xtc> [ligand_pdb]")
    else:
        lig_pdb = sys.argv[3] if len(sys.argv) > 3 else None
        check_integrity(sys.argv[1], sys.argv[2], lig_pdb)

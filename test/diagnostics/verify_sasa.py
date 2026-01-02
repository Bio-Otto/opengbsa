import mdtraj as md
import numpy as np

def verify_sasa():
    print("Verifying SASA (Solvent Accessible Surface Area)...")
    
    # Load trajectory
    traj_file = 'test/3if6_test/3if6.xtc'
    pdb_file = 'test/3if6_test/snapshot.pdb'
    
    print(f"Loading {traj_file} (Frame 0 only)...")
    # Optimize: Load only the first frame index 0
    traj = md.load_frame(traj_file, 0, top=pdb_file)
    frame = traj # It is already a single frame trajectory
    
    # Selections
    protein_sel = frame.topology.select("protein")
    ligand_sel = frame.topology.select("resname UNL")
    # FIX: Select only Protein + Ligand (exclude Water/Ions)
    complex_sel = frame.topology.select("protein or resname UNL")
    
    if len(ligand_sel) == 0:
        print("ERROR: Ligand (UNL) not found in topology!")
        return
        
    # Create sub-trajectories
    prot_traj = frame.atom_slice(protein_sel)
    lig_traj = frame.atom_slice(ligand_sel)
    comp_traj = frame.atom_slice(complex_sel)
    
    # Calculate SASA (Shrake-Rupley)
    # mode='residue' sums up for us? No, returns per-atom.
    # We want total SASA.
    
    sasa_prot_per_atom = md.shrake_rupley(prot_traj, mode='atom')
    sasa_lig_per_atom = md.shrake_rupley(lig_traj, mode='atom')
    sasa_comp_per_atom = md.shrake_rupley(comp_traj, mode='atom')
    
    total_sasa_prot = np.sum(sasa_prot_per_atom)
    total_sasa_lig = np.sum(sasa_lig_per_atom)
    total_sasa_comp = np.sum(sasa_comp_per_atom)
    
    # Convert nm^2 to Angstrom^2 (MDTraj uses nm^2, SASA is usually reported in A^2)
    # 1 nm^2 = 100 A^2
    sasa_p_ang = total_sasa_prot * 100
    sasa_l_ang = total_sasa_lig * 100
    sasa_c_ang = total_sasa_comp * 100
    
    print("-" * 30)
    print(f"Protein SASA: {sasa_p_ang:.2f} Å²")
    print(f"Ligand SASA:  {sasa_l_ang:.2f} Å²")  
    print(f"Complex SASA: {sasa_c_ang:.2f} Å²")
    print("-" * 30)
    
    sum_parts = sasa_p_ang + sasa_l_ang
    diff = sum_parts - sasa_c_ang
    
    print(f"Sum (Prot+Lig): {sum_parts:.2f} Å²")
    print(f"Buried Area:    {diff:.2f} Å²")
    
    if sasa_c_ang < sum_parts:
        print("✅ VERIFIED: Complex SASA < (Protein SASA + Ligand SASA)")
        print("   Physics holds true (Interface area is buried).")
    else:
        print("❌ FAILED: Complex SASA >= (Protein SASA + Ligand SASA)")
        print("   Something is wrong with the geometry or definition.")

if __name__ == "__main__":
    verify_sasa()

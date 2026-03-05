
import mdtraj as md
import numpy as np
import parmed as pmd

def calculate_sasa_energy(struct, probe_radius=1.4, surf_tension=0.0072):
    # Convert ParmEd to MDTraj
    # Save as PDB and reload is safest
    struct.save('temp_for_sasa.pdb', overwrite=True)
    traj = md.load('temp_for_sasa.pdb')
    
    # Calculate SASA (Shrake-Rupley)
    # mode='residue' gives per-residue, 'atom' gives per-atom
    sasa_per_atom = md.shrake_rupley(traj, probe_radius=probe_radius, mode='atom')
    
    # Sum total SASA (nm^2)
    sasa_total_nm2 = np.sum(sasa_per_atom)
    
    # Convert to Angstrom^2 if needed?
    # MDTraj returns nm^2.
    # OpenMM/Amber usually use Surface Tension in kcal/mol/A^2 or kJ/mol/nm^2
    # 0.0072 is typically kcal/mol/A^2.
    # 1 nm^2 = 100 A^2.
    
    sasa_total_A2 = sasa_total_nm2 * 100.0
    energy = sasa_total_A2 * surf_tension
    
    return sasa_total_A2, energy

print("Loading structures...")
# We need to calculate Delta SASA = Complex - Receptor - Ligand
# To correspond to the trajectory frame, we need the coordinates.
# Let's use the first frame from the generated PDBs if available, or just the input files.
# But input files (complex.pdb from data) might not match the frame analyzed.
# The previous analysis output `temp_from_traj.pdb` (Frame 0 usually).

try:
    complex_pdb = 'test/results/6t1h_1151_comp/analysis_final_fix/mmgbsa_analysis_20260208_185025/temp_from_traj.pdb'
    print(f"Using trajectory frame: {complex_pdb}")
    traj = md.load(complex_pdb)
except:
    print("Trajectory frame not found, using input pdb")
    complex_pdb = 'test/data/6t1h_1151_comp/complex.pdb'
    traj = md.load(complex_pdb)

# Load topology to identify receptor/ligand atoms
# Assuming 'LIG' resname
ligand_indices = traj.topology.select("resname LIG")
receptor_indices = traj.topology.select("not resname LIG")

print(f"Receptor atoms: {len(receptor_indices)}")
print(f"Ligand atoms: {len(ligand_indices)}")

# Create trajectories for components
traj_complex = traj
traj_receptor = traj.atom_slice(receptor_indices)
traj_ligand = traj.atom_slice(ligand_indices)

# Calculate SASA for each
# Using 0.0072 kcal/mol/A^2 (standard Amber/MMPBSA for LCPO?)
# Actually MMPBSA often uses 0.0072 or 0.0054 depending on method.
SURF_TENSION = 0.0072 

# Faster approx
POINTS = 50

print("Calculating SASA for COMPLEX...", flush=True)
sasa_c_nm2 = md.shrake_rupley(traj_complex, mode='atom', n_sphere_points=POINTS).sum()
print("Calculating SASA for RECEPTOR...", flush=True)
sasa_r_nm2 = md.shrake_rupley(traj_receptor, mode='atom', n_sphere_points=POINTS).sum()
print("Calculating SASA for LIGAND...", flush=True)
sasa_l_nm2 = md.shrake_rupley(traj_ligand, mode='atom', n_sphere_points=POINTS).sum()

sasa_c_A2 = sasa_c_nm2 * 100.0
sasa_r_A2 = sasa_r_nm2 * 100.0
sasa_l_A2 = sasa_l_nm2 * 100.0

e_c = sasa_c_A2 * SURF_TENSION
e_r = sasa_r_A2 * SURF_TENSION
e_l = sasa_l_A2 * SURF_TENSION

delta_sasa_A2 = sasa_c_A2 - sasa_r_A2 - sasa_l_A2
delta_e = e_c - e_r - e_l

print("-" * 50)
print(f"{'Component':<15} {'Area (A^2)':<15} {'Energy (kcal/mol)':<20}")
print("-" * 50)
print(f"{'Complex':<15} {sasa_c_A2:<15.2f} {e_c:<20.2f}")
print(f"{'Receptor':<15} {sasa_r_A2:<15.2f} {e_r:<20.2f}")
print(f"{'Ligand':<15} {sasa_l_A2:<15.2f} {e_l:<20.2f}")
print("-" * 50)
print(f"{'DELTA':<15} {delta_sasa_A2:<15.2f} {delta_e:<20.2f}")
print("-" * 50)
print(f"Using Surface Tension: {SURF_TENSION} kcal/mol/A^2")
print(f"Ref MMPBSA.py SA: -2.68 kcal/mol")

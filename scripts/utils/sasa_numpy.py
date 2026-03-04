
import numpy as np
import parmed as pmd
import sys

def generate_sphere_points(n):
    """
    Golden Section Spiral to generate N points on a sphere.
    """
    indices = np.arange(0, n, dtype=float) + 0.5
    phi = np.arccos(1 - 2*indices/n)
    theta = np.pi * (1 + 5**0.5) * indices
    x, y, z = np.cos(theta) * np.sin(phi), np.sin(theta) * np.sin(phi), np.cos(phi)
    return np.stack((x, y, z), axis=-1)

def calculate_sasa_numpy(coords, radii, probe_radius=1.4, n_points=96):
    """
    Calculate SASA using Shrake-Rupley algorithm (vectorized numpy).
    coords: (N, 3)
    radii: (N,)
    """
    n_atoms = len(coords)
    sphere_points = generate_sphere_points(n_points) # (M, 3)
    
    total_area = 0.0
    atom_areas = np.zeros(n_atoms)
    
    # Pre-calculate extended radii
    extended_radii = radii + probe_radius
    
    # Neighbor search (naive, O(N^2))
    # Optimize: Distance matrix slice?
    # For 4000 atoms, N^2 matrix is 16M floats ~ 128MB. Feasible.
    
    print("  Computing distance matrix...", flush=True)
    d2 = np.sum((coords[:, np.newaxis, :] - coords[np.newaxis, :, :])**2, axis=-1)
    
    print("  Calculating areas per atom...", flush=True)
    # Iterate atoms (doing all at once might be too heavy with broadcasting M points against N neighbors)
    for i in range(n_atoms):
        if i % 500 == 0: print(f"    Atom {i}/{n_atoms}", flush=True)
        
        r_i = extended_radii[i]
        pos_i = coords[i]
        
        # Find neighbors within cutoff (max possible radius + r_i)
        # Conservative cutoff: r_i + max(r_j)
        # Let's just use 10.0 A as safe cutoff for typical biomolecules
        cutoff_dist = r_i + np.max(extended_radii) 
        
        neighbor_indices = np.where(d2[i] < cutoff_dist**2)[0]
        neighbor_indices = neighbor_indices[neighbor_indices != i]
        
        if len(neighbor_indices) == 0:
            atom_areas[i] = 4 * np.pi * r_i**2
            continue
            
        # Points on sphere i
        points_i = pos_i + sphere_points * r_i # (M, 3)
        
        # Check against neighbors
        # For each point, is it inside any neighbor j?
        # dist(p, pos_j) < r_j
        
        # Broadcast: (M, 1, 3) - (1, K, 3) -> (M, K, 3)
        # This can be memory intensive if M and K are large.
        # Loop over batches of neighbors?
        # Or loop over neighbors?
        
        is_accessible = np.ones(n_points, dtype=bool)
        
        neighbors_pos = coords[neighbor_indices] # (K, 3)
        neighbors_r2 = extended_radii[neighbor_indices]**2 # (K,)
        
        # Check points against neighbors
        # Optimize: Check if point is inside ANY neighbor
        # (p - n)**2 < r_n**2
        
        for j in range(len(neighbor_indices)):
            nj_pos = neighbors_pos[j]
            nj_r2 = neighbors_r2[j]
            
            d2_points = np.sum((points_i - nj_pos)**2, axis=1)
            buried_mask = d2_points < nj_r2 - 1e-3 # Epsilon for stability
            is_accessible[buried_mask] = False
            
            if not np.any(is_accessible):
                break
                
        accessible_fraction = np.sum(is_accessible) / n_points
        area = 4 * np.pi * r_i**2 * accessible_fraction
        atom_areas[i] = area
        
    return np.sum(atom_areas)

prmtop_path = 'test/data/6t1h_1151_comp/complex.prmtop'
pdb_path = 'test/results/6t1h_1151_comp/analysis_final_fix/mmgbsa_analysis_20260208_185025/temp_from_traj.pdb'

print("Loading PDB for coordinates...", flush=True)
pdb = pmd.load_file(pdb_path)
coords = pdb.coordinates # (N, 3) 
if isinstance(coords, list): coords = coords[0]
if coords.shape[0] == 3 and coords.shape[1] > 3: coords = coords.T

print("Loading PRMTOP for radii...", flush=True)
structure = pmd.load_file(prmtop_path)
structure.coordinates = coords


radii = []
# ParmEd atoms have .rmin or .sigma?
# Amber prmtops usually have Rmin/2 or Sigma.
# atom.radii might be populated if radii set was used?
# Let's check atom.solvent_radius or simply deduce from sigma.
# Rmin = 2^(1/6) * sigma
# Radius = Rmin / 2
# But prmtop usually stores sigma/epsilon for LJ.
# For GB/SA, it stores specific PB/GB radii (atom.screen?).
# Let's try to find 'radius' attribute or 'solvent_radius'.
# If not, use LJ Rmin/2.

has_radii = False
for atom in structure.atoms:
    if hasattr(atom, 'solvent_radius'):
         # radii is usually intrinsic radius. solvent_radius might include probe?
         # No, typically just radius.
         radii.append(atom.solvent_radius)
         has_radii = True
    elif hasattr(atom, 'radius'):
         radii.append(atom.radius)
         has_radii = True
    else:
         # Fallback to LJ
         # eps, rmin = atom.epsilon, atom.rmin
         # rmin is usually 2*r_vdw
         radii.append(atom.rmin / 2.0)
         
radii = np.array(radii)
print(f"Radii loaded. Mean: {np.mean(radii):.2f}, Max: {np.max(radii):.2f}")


# Identification
ligand_mask = np.array([a.residue.name == 'LIG' for a in pdb.atoms])
receptor_mask = ~ligand_mask

c_coords = coords
c_radii = radii

r_coords = coords[receptor_mask]
r_radii = radii[receptor_mask]

l_coords = coords[ligand_mask]
l_radii = radii[ligand_mask]

SURF_TENSION = 0.0072
POINTS = 100

print(f"Calculating Complex SASA ({len(c_coords)} atoms)...", flush=True)
sasa_c = calculate_sasa_numpy(c_coords, c_radii, n_points=POINTS)

print(f"Calculating Receptor SASA ({len(r_coords)} atoms)...", flush=True)
sasa_r = calculate_sasa_numpy(r_coords, r_radii, n_points=POINTS)

print(f"Calculating Ligand SASA ({len(l_coords)} atoms)...", flush=True)
sasa_l = calculate_sasa_numpy(l_coords, l_radii, n_points=POINTS)

delta_sasa = sasa_c - sasa_r - sasa_l
delta_e = delta_sasa * SURF_TENSION

print("-" * 50)
print(f"{'Component':<15} {'Area (A^2)':<15} {'Energy (kcal/mol)':<20}")
print("-" * 50)
print(f"{'Complex':<15} {sasa_c:<15.2f} {sasa_c * SURF_TENSION:<20.2f}")
print(f"{'Receptor':<15} {sasa_r:<15.2f} {sasa_r * SURF_TENSION:<20.2f}")
print(f"{'Ligand':<15} {sasa_l:<15.2f} {sasa_l * SURF_TENSION:<20.2f}")
print("-" * 50)
print(f"{'DELTA':<15} {delta_sasa:<15.2f} {delta_e:<20.2f}")
print("-" * 50)
print(f"Ref MMPBSA.py SA: -2.68 kcal/mol")

#!/usr/bin/env python3
"""
Check if the atom order in ligand2.sdf matches the UNL residue in complex2_fixed.pdb.
"""

from openff.toolkit.topology import Molecule
from openmm import app

sdf_file = 'test/ligand2.sdf'
pdb_file = 'test/complex2_fixed.pdb'

print(f"Loading SDF {sdf_file}...")
mol_sdf = Molecule.from_file(sdf_file)
print(f"SDF Atoms: {mol_sdf.n_atoms}")

print(f"Loading PDB {pdb_file}...")
pdb = app.PDBFile(pdb_file)
# Extract UNL atoms
unl_atoms = [a for a in pdb.topology.atoms() if a.residue.name == 'UNL']
print(f"PDB UNL Atoms: {len(unl_atoms)}")

if mol_sdf.n_atoms != len(unl_atoms):
    print("ERROR: Atom count mismatch!")
    exit(1)

# To check order, we need to create an OpenFF Molecule from the PDB atoms
# But PDB lacks bond orders.
# However, we can infer partial isomorphism or use RDKit.

# Better: Try to assign the PDB positions to the SDF molecule based on simple order
# and see if it makes geometric sense (bond lengths).
# If atoms are scrambled, we'll see huge bond lengths.

import numpy as np
from openmm import unit

positions = pdb.positions
# Extract UNL positions in PDB order
unl_positions = []
for a in unl_atoms:
    pos = positions[a.index].value_in_unit(unit.angstrom)
    unl_positions.append(pos)
unl_positions = np.array(unl_positions)

# SDF Connectivity
# Iterate over bonds in SDF and calculate distance in PDB positions
print("\nChecking bond lengths assuming 1:1 mapping...")
max_bond_len = 0.0
deviations = 0
for bond in mol_sdf.bonds:
    idx1 = bond.atom1_index
    idx2 = bond.atom2_index
    
    pos1 = unl_positions[idx1]
    pos2 = unl_positions[idx2]
    
    dist = np.linalg.norm(pos1 - pos2)
    if dist > 2.0: # 2 Angstrom threshold for bonded
        print(f"  Bond {idx1}-{idx2}: Distance {dist:.2f} A (Expected ~1.5 A) -> MISMATCH")
        deviations += 1
    max_bond_len = max(max_bond_len, dist)

print(f"Max bond length with 1:1 mapping: {max_bond_len:.2f} A")
if deviations > 0:
    print(f"Found {deviations} bonds with bad lengths. CONCLUSION: ATOM ORDER MISMATCH.")
else:
    print("Bond lengths look reasonable. Order might be correct.")

if deviations > 0:
    print("\nAttempting to find correct mapping...")
    # Map PDB topology (connectivity) to SDF topology?
    # PDB has bonds (detected by OpenMM based on distance or standard residues).
    # Since UNL is not standard, OpenMM might compute bonds based on distance.
    # Let's verify PDB bonds.
    
    unl_pdb_bonds = []
    # Map PDB atom index to 0..N local index
    pdb_idx_map = {a.index: i for i, a in enumerate(unl_atoms)}
    
    for bond in pdb.topology.bonds():
        if bond.atom1.residue.name == 'UNL' and bond.atom2.residue.name == 'UNL':
            i1 = pdb_idx_map[bond.atom1.index]
            i2 = pdb_idx_map[bond.atom2.index]
            unl_pdb_bonds.append(sorted((i1, i2)))
    
    print(f"PDB identified {len(unl_pdb_bonds)} bonds within UNL based on distance.")
    
    # We can try to remap using OpenFF isomorphic match if we can construct a graph from PDB.
    # But simpler: we know the PDB positions are valid (user says VMD looks good).
    # We want to reorder PDB atoms to match SDF order.
    
    # Use OpenFF to match
    # Create mol from PDB (using connectivity guess)
    try:
        mol_pdb = Molecule.from_topology(pdb.topology) # Hard without atoms info
        # Actually Molecule.from_openmm works if we allow partial info
        pass 
    except:
        pass
        
    # Isomorphism is hard without element matching.
    # Let's assume Elements are correct but shuffled? 
    # Or just check if we can remap based on Graph Isomorphism (NetworkX).
    
    # Actually, RDKit is best for this.
    try:
        from rdkit import Chem
        # Create RDKit mol from PDB block?
        # or just atoms/bonds
        # We need to construct RDKit mol from SDF and PDB and match.
        pass
    except ImportError:
        print("RDKit not installed?")


#!/usr/bin/env python3
"""
Check for steric clashes in complex2_fixed.pdb.
Clash defined as distance < sum of radii * 0.7?
Or simply < 1.0 A for Heavy-Heavy.
"""

from openmm import app, unit
import numpy as np
import itertools

pdb_file = 'test/complex2_fixed.pdb'
print(f"Loading {pdb_file}...")
pdb = app.PDBFile(pdb_file)

atoms = list(pdb.topology.atoms())
positions = pdb.positions.value_in_unit(unit.angstrom)

ligand_atoms = [a for a in atoms if a.residue.name == 'UNL']
protein_atoms = [a for a in atoms if a.residue.name != 'UNL' and a.residue.name != 'HOH']

print(f"Ligand atoms: {len(ligand_atoms)}")
print(f"Protein atoms: {len(protein_atoms)}")

# Naive O(N*M) check? 8000*138 ~ 1 million. Fast enough.

print("Checking Heavy-Heavy clashes (dist < 2.0 A)...")
clashes = 0
min_dist = 999.0
clash_pair = None

for la in ligand_atoms:
    if la.element.name == 'hydrogen': continue
    l_pos = positions[la.index]
    
    for pa in protein_atoms:
        if pa.element.name == 'hydrogen': continue
        p_pos = positions[pa.index]
        
        dist = np.linalg.norm(l_pos - p_pos)
        min_dist = min(min_dist, dist)
        
        if dist < 2.0:
            print(f"CLASH! {la.name} (Lig) - {pa.name} ({pa.residue.name} {pa.residue.id}): {dist:.2f} A")
            clashes += 1
            if clashes > 10:
                print("... truncated ...")
                break
    if clashes > 10: break

print(f"Min Heavy-Heavy Dist: {min_dist:.2f} A")

print("\nChecking H-Heavy or H-H clashes (dist < 1.5 A assuming H)...")
# Check scans involving H
h_clashes = 0
min_h_dist = 999.0

# Optimize: Put protein atoms in grid? No need for 1M check.
# Just run it.

for la in ligand_atoms:
    l_pos = positions[la.index]
    l_is_h = (la.element.name == 'hydrogen')
    
    for pa in protein_atoms:
        p_pos = positions[pa.index]
        p_is_h = (pa.element.name == 'hydrogen')
        
        if not l_is_h and not p_is_h: continue # Already checked
        
        threshold = 1.5 # Strict
        dist = np.linalg.norm(l_pos - p_pos)
        min_h_dist = min(min_h_dist, dist)
        
        if dist < threshold:
             # print(f"H-CLASH! {la.name} - {pa.name}: {dist:.2f} A")
             h_clashes += 1

print(f"Total H-involved clashes (< 1.5 A): {h_clashes}")
print(f"Min H-involved Dist: {min_h_dist:.2f} A")

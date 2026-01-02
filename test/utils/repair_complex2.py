#!/usr/bin/env python3
"""
Repair complex2 PDB using Modeller on Protein component.
1. Extract Protein.
2. Fix Protein (addHydrogens).
3. Extract Ligand.
4. Merge.
"""

from openmm import app
import openmm
from openmm import unit

input_pdb = 'test/complex2.pdb'
output_pdb = 'test/complex2_fixed.pdb'
ligand_resname = 'UNL'

print(f"Loading {input_pdb}...")
pdb = app.PDBFile(input_pdb)

# Separation using Modeller
modeller = app.Modeller(pdb.topology, pdb.positions)

# Identify ligand residues to remove for protein repair
# We want to keep ONLY protein first
ligand_residues = [r for r in modeller.topology.residues() if r.name == ligand_resname]
modeller.delete(ligand_residues)

print(f"Protein extraction: {modeller.topology.getNumAtoms()} atoms.")

# Repair Protein (add hydrogens, fix terminals)
print("Adding hydrogens to protein...")
forcefield = app.ForceField('amber14-all.xml', 'amber14/tip3p.xml')
modeller.addHydrogens(forcefield)

protein_topology = modeller.topology
protein_positions = modeller.positions

print(f"Repaired protein: {protein_topology.getNumAtoms()} atoms.")

# Now we need to get the ligand back.
# Reload original PDB to get ligand
pdb_full = app.PDBFile(input_pdb)
modeller_full = app.Modeller(pdb_full.topology, pdb_full.positions)
# Delete everything NOT ligand
protein_residues = [r for r in modeller_full.topology.residues() if r.name != ligand_resname]
modeller_full.delete(protein_residues)

ligand_topology = modeller_full.topology
ligand_positions = modeller_full.positions
print(f"Ligand extraction: {ligand_topology.getNumAtoms()} atoms.")

# Merge
# Modeller.add invokes data copy.
# We treat the repaired protein as base.
# We create a new Modeller with repaired protein, then add ligand.
# Wait, Modeller.add(hessian/topology...) - NO.
# Modeller.add(other_topology, other_positions)
# Let's check Modeller source or docs.
# Modeller.add(topology, positions) joins them.

print("Merging repaired protein and ligand...")
# Create fresh modeller with protein
final_modeller = app.Modeller(protein_topology, protein_positions)
final_modeller.add(ligand_topology, ligand_positions)

print(f"Final system: {final_modeller.topology.getNumAtoms()} atoms.")

print(f"Saving to {output_pdb}...")
with open(output_pdb, 'w') as f:
    app.PDBFile.writeFile(final_modeller.topology, final_modeller.positions, f, keepIds=True)

print("Done!")

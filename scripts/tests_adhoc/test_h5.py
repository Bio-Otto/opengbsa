import openmm.app as app
import openmm.unit as unit

pdb = app.PDBFile('test/data/3if6_dimer_test/3if6_sanitized.pdb')

modeller = app.Modeller(pdb.topology, pdb.positions)
to_delete_res = [r for r in modeller.topology.residues() if r.name in ['HOH', 'WAT', 'TIP3', 'SOL', 'NA', 'CL', 'K', 'MG', 'ZN', 'CA', 'LIG', 'UNL']]
modeller.delete(to_delete_res)

# Find cross-chain bonds
cross_chain_bonds = []
for a1, a2 in modeller.topology.bonds():
    if a1.residue.chain != a2.residue.chain:
        cross_chain_bonds.append((a1, a2))
        
print(f"Found {len(cross_chain_bonds)} cross-chain bonds:")
for a1, a2 in cross_chain_bonds:
    print(f"  {a1.residue.name}{a1.residue.index}.{a1.name} (Chain {a1.residue.chain.index}) -- {a2.residue.name}{a2.residue.index}.{a2.name} (Chain {a2.residue.chain.index})")

# Remove them? Modeller doesn't have a delete_bond. We have to reconstruct the topology or delete the atoms.
# Wait, let's see why PDBFile added a bond between chains.

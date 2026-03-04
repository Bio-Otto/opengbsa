import openmm.app as app
import openmm.unit as unit

pdb = app.PDBFile('test/data/3if6_dimer_test/3if6_sanitized.pdb')

modeller = app.Modeller(pdb.topology, pdb.positions)
to_delete_res = [r for r in modeller.topology.residues() if r.name in ['HOH', 'WAT', 'TIP3', 'SOL', 'NA', 'CL', 'K', 'MG', 'ZN', 'CA', 'LIG', 'UNL']]
modeller.delete(to_delete_res)

# Let's inspect ASN 241
for r in modeller.topology.residues():
    if r.index == 241:
        print(f"Residue {r.name} {r.index}:")
        for a in r.atoms():
            print(f"  Atom {a.name} ({a.element}) bonded to:", [b.name for (a1, a2) in modeller.topology.bonds() if a1==a or a2==a for b in (a1, a2) if b!=a])


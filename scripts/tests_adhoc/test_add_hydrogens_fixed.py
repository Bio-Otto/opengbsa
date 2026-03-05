from openmm import app

pdb = app.PDBFile('test/data/3if6_dimer_test/3if6_fixed.pdb')
print("PDB Topology loaded from cleaned file.")

forcefield = app.ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')

print("Modeller...")
modeller = app.Modeller(pdb.topology, pdb.positions)
modeller.addHydrogens(forcefield)

print("Add hydrogens successful!")

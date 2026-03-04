import parmed as pmd
from openmm import app

print("Loading with ParmEd...")
struct = pmd.load_file('test/data/3if6_dimer_test/3if6_fixed.pdb')
omm_top = struct.topology
omm_pos = struct.positions

print("Creating Modeller from ParmEd's topology...")
modeller = app.Modeller(omm_top, omm_pos)

forcefield = app.ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
print("Adding Hydrogens...")
modeller.addHydrogens(forcefield)

print("Add hydrogens successful!")

import openmm.app as app
import openmm.unit as unit
from openmmforcefields.generators import SystemGenerator

pdb = app.PDBFile('test/data/3if6_dimer_test/3if6_sanitized.pdb')
print("PDB LOADED. Atoms:", pdb.topology.getNumAtoms())

modeller = app.Modeller(pdb.topology, pdb.positions)
to_delete_res = [r for r in modeller.topology.residues() if r.name in ['HOH', 'WAT', 'TIP3', 'SOL', 'NA', 'CL', 'K', 'MG', 'ZN', 'CA', 'LIG', 'UNL']]
modeller.delete(to_delete_res)

to_delete_h = [atom for atom in modeller.topology.atoms() if atom.element is not None and atom.element.symbol == 'H']
if to_delete_h:
    modeller.delete(to_delete_h)
    print(f"Stripped {len(to_delete_h)} hydrogens")

sysgen = SystemGenerator(forcefields=['amber/ff14SB.xml', 'amber/tip3p_standard.xml'])

try:
    modeller.addHydrogens(sysgen.forcefield)
    print("Hydrogens added successfully!")
except Exception as e:
    print("addHydrogens ERROR:", type(e).__name__, e)

try:
    sys = sysgen.create_system(modeller.topology)
    print("System created successfully!")
except Exception as e:
    print("create_system ERROR:", type(e).__name__, e)

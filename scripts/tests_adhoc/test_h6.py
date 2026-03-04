import openmm.app as app
import openmm.unit as unit

pdb = app.PDBFile('test/data/3if6_dimer_test/3if6_sanitized.pdb')

modeller = app.Modeller(pdb.topology, pdb.positions)
to_delete_res = [r for r in modeller.topology.residues() if r.name in ['HOH', 'WAT', 'TIP3', 'SOL', 'NA', 'CL', 'K', 'MG', 'ZN', 'CA', 'LIG', 'UNL']]
modeller.delete(to_delete_res)

# We want to delete ALL hydrogens ONLY from the N-terminal residues (VAL 0 and ILE 242)
# and let addHydrogens strictly protonate them according to NVAL/NILE.
termini_indices = [0, 242]
to_delete_h = []
for res in modeller.topology.residues():
    if res.index in termini_indices:
        for a in res.atoms():
            if a.element is not None and a.element.symbol == 'H':
                to_delete_h.append(a)

if to_delete_h:
    modeller.delete(to_delete_h)
    print(f"Deleted {len(to_delete_h)} hydrogens from termini")

forcefield = app.ForceField('amber/ff14SB.xml', 'amber/tip3p_standard.xml')

# Try addHydrogens with default
try:
    modeller.addHydrogens(forcefield)
    print("Hydrogens added to termini successfully!")
except Exception as e:
    print("addHydrogens ERROR:", type(e).__name__, e)

try:
    sys = forcefield.createSystem(modeller.topology, nonbondedMethod=app.NoCutoff, constraints=app.HBonds)
    print("System created successfully!")
except Exception as e:
    print("create_system ERROR:", type(e).__name__, e)


import openmm.app as app
import openmm.unit as unit

pdb = app.PDBFile('test/data/3if6_dimer_test/3if6_sanitized.pdb')
print("PDB LOADED")

modeller = app.Modeller(pdb.topology, pdb.positions)
to_delete_res = [r for r in modeller.topology.residues() if r.name in ['HOH', 'WAT', 'TIP3', 'SOL', 'NA', 'CL', 'K', 'MG', 'ZN', 'CA', 'LIG', 'UNL']]
modeller.delete(to_delete_res)
print("Removed water/ions.")

forcefield = app.ForceField('amber/ff14SB.xml', 'amber/tip3p_standard.xml')

# Try createSystem with ignoreExternalBonds=True without ANY addHydrogens or stripping
try:
    sys = forcefield.createSystem(
        modeller.topology, 
        nonbondedMethod=app.NoCutoff, 
        constraints=app.HBonds, 
        ignoreExternalBonds=True
    )
    print("System created successfully with HBonds constraints!")
except Exception as e:
    print("create_system (HBonds) ERROR:", type(e).__name__, e)

try:
    sys2 = forcefield.createSystem(
        modeller.topology, 
        nonbondedMethod=app.NoCutoff, 
        constraints=None, 
        ignoreExternalBonds=True
    )
    print("System created successfully with NO constraints!")
except Exception as e:
    print("create_system (No constraints) ERROR:", type(e).__name__, e)

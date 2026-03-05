import simtk.openmm.app as app
import simtk.openmm as openmm
import simtk.unit as unit

print("Loading 3if6.pdb...")
pdb = app.PDBFile('test/data/3if6_dimer_test/3if6.pdb')
modeller = app.Modeller(pdb.topology, pdb.positions)

ff = app.ForceField('amber/ff14SB.xml', 'amber/tip3p_standard.xml')

print("Generating system directly WITHOUT addHydrogens...")
try:
    system = ff.createSystem(modeller.topology)
    print("Success without addHydrogens!")
except Exception as e:
    print("Failed without addHydrogens:", e)

print("Attempting to addHydrogens()...")
try:
    modeller.addHydrogens(ff)
    print("addHydrogens() successful. Num atoms:", modeller.topology.getNumAtoms())
    print("Generating system...")
    system_fixed = ff.createSystem(modeller.topology)
    print("Success!")
except Exception as e:
    print("Error during addHydrogens or createSystem:", e)

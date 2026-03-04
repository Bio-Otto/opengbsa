
import sys
import os
from openmm import app, unit, openmm

def get_energy(system, positions):
    integrator = openmm.VerletIntegrator(1.0*unit.femtoseconds)
    context = openmm.Context(system, integrator)
    context.setPositions(positions)
    # Minimize slightly
    openmm.LocalEnergyMinimizer.minimize(context, maxIterations=50)
    state = context.getState(getEnergy=True)
    return state.getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole)

def run_test():
    pdb_path = 'test/data/6xj3_pdb_test/6xj3.pdb'
    if not os.path.exists(pdb_path):
        print("Missing pdb file.")
        return

    print("Loading PDB...")
    pdb = app.PDBFile(pdb_path)
    modeller = app.Modeller(pdb.topology, pdb.positions)
    
    # Keep only protein (standard residues)
    standard_res = {'ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL',
                    'HID','HIE','HIP','CYX','ASH','GLH','LYN'}
    
    to_delete = [r for r in modeller.topology.residues() if r.name not in standard_res]
    if to_delete:
        print(f"Removing {len(to_delete)} non-standard residues (ligands/waters/ions).")
        modeller.delete(to_delete)
        
    print(f"Remaining atoms: {modeller.topology.getNumAtoms()}")
    
    
    # OBC1
    print("\n--- Testing OBC1 (PDB Mode + XML) ---")
    ff1 = app.ForceField('amber14-all.xml', 'amber14/tip3p.xml', 'implicit/obc1.xml')
    sys1 = ff1.createSystem(
        modeller.topology,
        nonbondedMethod=app.NoCutoff,
        constraints=app.HBonds
    )
    # Emulate refine (SA=0)
    for f in sys1.getForces():
        if isinstance(f, openmm.GBSAOBCForce):
            f.setSurfaceAreaEnergy(0.0)

    e1 = get_energy(sys1, modeller.positions)
    print(f"OBC1 Energy: {e1:.4f} kcal/mol")
    
    # OBC2
    print("\n--- Testing OBC2 (PDB Mode + XML) ---")
    ff2 = app.ForceField('amber14-all.xml', 'amber14/tip3p.xml', 'implicit/obc2.xml')
    sys2 = ff2.createSystem(
        modeller.topology,
        nonbondedMethod=app.NoCutoff,
        constraints=app.HBonds
    )
    # Emulate refine (SA=0)
    for f in sys2.getForces():
         if isinstance(f, openmm.GBSAOBCForce):
            f.setSurfaceAreaEnergy(0.0)

    e2 = get_energy(sys2, modeller.positions)
    print(f"OBC2 Energy: {e2:.4f} kcal/mol")

    
    diff = abs(e1 - e2)
    print(f"\nEnergy Difference: {diff:.4f} kcal/mol")
    
    if diff > 1.0: # Expect large difference for whole protein
        print("SUCCESS: OBC1 and OBC2 produce different energies via ForceField!")
    else:
        print("FAILURE: Energies are identical.")

if __name__ == "__main__":
    run_test()

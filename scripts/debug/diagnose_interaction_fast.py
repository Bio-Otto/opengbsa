
import parmed as pmd
from openmm import app, openmm, unit
import mdtraj as md
import numpy as np
import time

def calculate_interaction_fast():
    print("Starting Fast Interaction Diagnostic (NoCutoff)...")
    prmtop = "test/data/6t1h_6466_comp/complex.prmtop"
    xtc = "test/data/6t1h_6466_comp/md_complex_prod.xtc"
    solv_top = "test/data/6t1h_6466_comp/replica_1/gro/complex_solv_ions.gro"
    
    struct = pmd.load_file(prmtop)
    traj = md.load_frame(xtc, 0, top=solv_top)
    
    # Simple Slicing
    selection = np.arange(len(struct.atoms))
    frame = traj.atom_slice(selection)
    positions = frame.xyz[0] * unit.nanometer
    
    # System
    system = struct.createSystem(nonbondedMethod=app.NoCutoff)
    
    # Separate Forces
    nb = [f for f in system.getForces() if isinstance(f, openmm.NonbondedForce)][0]
    nb.setForceGroup(0)
    
    # Masks
    ligand_resname = "LIG"
    ligand_mask = [a.idx for a in struct.atoms if a.residue.name == ligand_resname]
    receptor_mask = [a.idx for a in struct.atoms if a.residue.name != ligand_resname]
    
    integrator = openmm.VerletIntegrator(1.0*unit.femtoseconds)
    ctx = openmm.Context(system, integrator)
    ctx.setPositions(positions)
    
    # 1. Complex Energy
    e_complex = ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole)
    print(f"Complex Potential: {e_complex:.2f}")
    
    # 2. Receptor Energy (Zero Ligand)
    # Save original parameters for ligand
    lig_params = {}
    for i in ligand_mask:
        lig_params[i] = nb.getParticleParameters(i)
        nb.setParticleParameters(i, 0.0, 0.1, 0.0)
        
    ctx.reinitialize(preserveState=True)
    e_receptor = ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole)
    print(f"Receptor Potential: {e_receptor:.2f}")
    
    # Restore Ligand
    for i, params in lig_params.items():
        nb.setParticleParameters(i, *params)
        
    # 3. Ligand Energy (Zero Receptor)
    for i in receptor_mask:
        nb.setParticleParameters(i, 0.0, 0.1, 0.0)
        
    ctx.reinitialize(preserveState=True)
    e_ligand = ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole)
    print(f"Ligand Potential: {e_ligand:.2f}")
    
    delta = e_complex - e_receptor - e_ligand
    print(f"Interaction Delta (Zeroing Method): {delta:.4f} kcal/mol")

if __name__ == "__main__":
    calculate_interaction_fast()

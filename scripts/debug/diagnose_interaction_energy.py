
import parmed as pmd
from openmm import app, openmm, unit
import mdtraj as md
import numpy as np

def calculate_interaction_energy():
    print("Loading structure...")
    prmtop = "test/data/6t1h_6466_comp/complex.prmtop"
    xtc = "test/data/6t1h_6466_comp/md_complex_prod.xtc"
    solv_top = "test/data/6t1h_6466_comp/replica_1/gro/complex_solv_ions.gro"
    
    struct = pmd.load_file(prmtop)
    print("Loading trajectory with solvated topology...")
    traj = md.load(xtc, top=solv_top)
    
    # Slice to match Dry Prmtop (Protein + Ligand)
    # Assuming standard selection "protein or resname LIG"
    print("Slicing trajectory to match Dry Complex...")
    selection = traj.topology.select("protein or resname LIG")
    if len(selection) != len(struct.atoms):
        print(f"Warning: Selection count {len(selection)} != Prmtop count {len(struct.atoms)}")
        # Try simplistic first N atoms if counts match that way (common in GROMACS -> Amber conversion)
        if len(selection) < len(struct.atoms):
             print("Selection too small, checking atom counts...")
        
    frame = traj[0].atom_slice(selection)
    positions = frame.xyz[0] * unit.nanometer
    
    # Create System (Vacuum, NoCutoff)
    system = struct.createSystem(nonbondedMethod=app.NoCutoff)
    
    # Split NonbondedForce to separate VDW/ELE and 1-4s
    for f in system.getForces():
        if isinstance(f, openmm.NonbondedForce):
            f.setForceGroup(0) # Standard NB
            
    # Calculate Total Energy
    integrator = openmm.VerletIntegrator(1.0*unit.femtoseconds)
    ctx = openmm.Context(system, integrator)
    ctx.setPositions(positions)
    
    state = ctx.getState(getEnergy=True)
    potential_e = state.getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole)
    print(f"Total Potential Energy (Complex): {potential_e:.4f} kcal/mol")
    
    # Define Masks
    ligand_resname = "LIG"
    ligand_mask = [a.index for a in struct.atoms if a.residue.name == ligand_resname]
    receptor_mask = [a.index for a in struct.atoms if a.residue.name != ligand_resname]
    
    print(f"Ligand Atoms: {len(ligand_mask)}")
    print(f"Receptor Atoms: {len(receptor_mask)}")
    
    # Calculate Interaction Energy Manually (Pairwise)
    # E_inter = Sum(q1*q2/r) + Sum(LJs)
    # This ignores 1-4s (which are Intra).
    # If this matches OpenGBSA Delta, then Delta is purely Inter.
    
    nb = [f for f in system.getForces() if isinstance(f, openmm.NonbondedForce)][0]
    
    elec_inter = 0.0
    vdw_inter = 0.0
    
    charges = []
    sigmas = []
    epsilons = []
    
    for i in range(nb.getNumParticles()):
        chg, sig, eps = nb.getParticleParameters(i)
        charges.append(chg.value_in_unit(unit.elementary_charge))
        sigmas.append(sig.value_in_unit(unit.nanometer))
        epsilons.append(eps.value_in_unit(unit.kilojoule_per_mole))
        
    pos = frame.xyz[0] # nm
    
    print("Calculating Pairwise Interaction (This is slow, analyzing subset or using OpenMM groups)...")
    
    # Use OpenMM Groups to separate Receptor and Ligand
    # This is faster and uses the exact ForceField logic
    
    for i in ligand_mask:
        nb.setParticleParameters(i, 0.0, 0.1, 0.0) # Zero out ligand
        
    ctx.reinitialize(preserveState=True)
    e_receptor = ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole)
    
    # Reset
    for i in range(nb.getNumParticles()):
        chg = charges[i]
        sig = sigmas[i] * unit.nanometer
        eps = epsilons[i] * unit.kilojoule_per_mole
        nb.setParticleParameters(i, chg, sig, eps)
        
    for i in receptor_mask:
        nb.setParticleParameters(i, 0.0, 0.1, 0.0) # Zero out receptor
        
    ctx.reinitialize(preserveState=True)
    e_ligand = ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole)
    
    delta = potential_e - e_receptor - e_ligand
    print(f"\nDecomposition (Total Potential):")
    print(f"  Complex:  {potential_e:.4f}")
    print(f"  Receptor: {e_receptor:.4f}")
    print(f"  Ligand:   {e_ligand:.4f}")
    print(f"  Delta:    {delta:.4f} kcal/mol")
    
    print("\nResult Analysis:")
    if abs(delta + 62) < 10:
        print("  Delta matches OpenGBSA result (~ -62).")
    else:
        print(f"  Delta ({delta}) differs from OpenGBSA (-62).")
        
if __name__ == "__main__":
    calculate_interaction_energy()

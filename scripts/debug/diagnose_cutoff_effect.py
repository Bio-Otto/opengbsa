
import parmed as pmd
from openmm import app, openmm, unit
import mdtraj as md
import numpy as np

def calculate_interaction_energy_cutoff(cutoff_dist):
    print(f"\n--- Testing Cutoff: {cutoff_dist} A ---")
    prmtop = "test/data/6t1h_6466_comp/complex.prmtop"
    xtc = "test/data/6t1h_6466_comp/md_complex_prod.xtc"
    solv_top = "test/data/6t1h_6466_comp/replica_1/gro/complex_solv_ions.gro"
    
    struct = pmd.load_file(prmtop)
    traj = md.load_frame(xtc, 0, top=solv_top)
    
    # Selection logic
    # We want to match the atoms in 'struct' (Protein + Ligand)
    # Usually they are the first N atoms in the solvated topology
    N_atoms = len(struct.atoms)
    
    # 1. Try selecting by index (0 to N-1)
    selection = np.arange(N_atoms)
    
    print(f"Checking consistency for first {N_atoms} atoms...")
    top = traj.topology
    # Check a few random atoms to verify alignment
    indices_to_check = [0, N_atoms//2, N_atoms-1]
    
    mismatch = False
    for i in indices_to_check:
        s_res = struct.atoms[i].residue.name
        t_res = top.atom(i).residue.name
        if s_res != t_res:
            # Map common differences (HOH vs WAT, etc. but here we expect PROTEIN)
            print(f"  Mismatch at {i}: Prmtop={s_res} vs Traj={t_res}")
            mismatch = True
            
    if mismatch:
        print("  Index mismatch detected! Trying to find matching atoms via name/residue...")
        # This is complex. Standard GROMACS to Amber often preserves order of Solute.
        # If mismatch, maybe Solute is not at start?
        pass # Proceed anyway to see result, but warn.

    frame = traj.atom_slice(selection)
    positions = frame.xyz[0] * unit.nanometer
    
    # Create System
    if cutoff_dist >= 999.0:
        method = app.NoCutoff
        cutoff = 999.0 * unit.angstroms
        print("Method: NoCutoff")
    else:
        method = app.CutoffNonPeriodic
        cutoff = cutoff_dist * unit.angstroms
        print(f"Method: CutoffNonPeriodic ({cutoff})")
        
    system = struct.createSystem(nonbondedMethod=method, nonbondedCutoff=cutoff)
    
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
    
    # Calculate Energies
    e_complex = ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole)
    
    # Zero Ligand
    original_params = []
    for i in ligand_mask:
        p = nb.getParticleParameters(i)
        original_params.append((i, p))
        nb.setParticleParameters(i, 0.0, 0.1, 0.0) 
    
    ctx.reinitialize(preserveState=True)
    e_receptor = ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole)
    
    # Restore Ligand
    for i, p in original_params:
        nb.setParticleParameters(i, *p)
        
    # Zero Receptor
    for i in receptor_mask:
        nb.setParticleParameters(i, 0.0, 0.1, 0.0)
        
    ctx.reinitialize(preserveState=True)
    e_ligand = ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole)
    
    delta = e_complex - e_receptor - e_ligand
    print(f"  Complex E: {e_complex:.2f}")
    print(f"  Receptor E: {e_receptor:.2f}")
    print(f"  Ligand E: {e_ligand:.2f}")
    print(f"  Interaction Delta: {delta:.4f} kcal/mol")
    
    return delta

if __name__ == "__main__":
    results = {}
    cutoffs = [999.0, 20.0, 12.0, 9.0]
    for c in cutoffs:
        d = calculate_interaction_energy_cutoff(c)
        results[c] = d
    
    print("\nSummary:")
    for c, d in results.items():
        print(f"Cutoff {c:>5} A: {d:.4f} kcal/mol")

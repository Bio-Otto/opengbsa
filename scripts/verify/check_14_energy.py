
import parmed as pmd
from openmm import app, openmm, unit
import mdtraj as md

def calculate_14_energy():
    print("Checking 1-4 Interaction Energy Magnitude...")
    prmtop = "test/data/6t1h_6466_comp/complex.prmtop"
    struct = pmd.load_file(prmtop)
    
    # Create System with NO Nonbonded, ONLY 1-4s (if possible)
    # Actually, create standard system, but zero out non-1-4s?
    # Easier: Iterate exceptions in NonbondedForce.
    
    system = struct.createSystem(nonbondedMethod=app.NoCutoff)
    nb = [f for f in system.getForces() if isinstance(f, openmm.NonbondedForce)][0]
    
    # We want to sum E_14 = q1q2/r + epsilon*((sigma/r)^12 - (sigma/r)^6)
    # But we need positions.
    # Load Frame 0 from XTC
    xtc = "test/data/6t1h_6466_comp/md_complex_prod.xtc"
    try:
        traj = md.load(xtc, top=prmtop)
        frame = traj[0]
        positions = frame.xyz[0] * unit.nanometer
    except Exception:
        print("Using solvated load logic...")
        solv_top = "test/data/6t1h_6466_comp/replica_1/gro/complex_solv_ions.gro"
        traj = md.load(xtc, top=solv_top)
        selection = traj.topology.select("protein or resname LIG")
        frame = traj[0].atom_slice(selection)
        positions = frame.xyz[0] * unit.nanometer

    # Calculate 1-4 Energy
    # Method: 
    # 1. Zero out all charges and epsilons for normal particles.
    # 2. Keep exceptions (which define 1-4s).
    # 3. Actually exceptions usually Scale relative to normal. 
    #    If normal is zero, 1-4 might be zeroed too depending on implementation?
    #    No, exceptions are explicit parameters in OpenMM.
    
    # Strategy: 
    # Set all main particle parameters to 0. Use exception parameters as is.
    # Verify if exception indices correspond to 1-4s.
    
    for i in range(nb.getNumParticles()):
        nb.setParticleParameters(i, 0.0, 0.1, 0.0)
        
    # Now E_nonbonded = Sum(Exceptions)
    # Note: Exclusions (1-2, 1-3) have chargeProd=0, epsilon=0.
    # 1-4s have non-zero.
    
    integrator = openmm.VerletIntegrator(1.0*unit.femtoseconds)
    ctx = openmm.Context(system, integrator)
    ctx.setPositions(positions)
    
    e_14 = ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole)
    print(f"Total 1-4 Energy (ELE+VDW): {e_14:.4f} kcal/mol")
    
    # Separate VDW and ELE 1-4?
    # Can't easily without rebuilding system twice.
    # But total magnitude gives us a clue.
    
    print("\nIf this value is large (e.g. > 100 kcal/mol), then 1-4 cancellation is critical.")
    print("If Complex 1-4s != Receptor + Ligand 1-4s exactly, artifacts appear.")

if __name__ == "__main__":
    calculate_14_energy()

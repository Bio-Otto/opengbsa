
import openmm
import openmm.app as app
import openmm.unit as unit
import numpy as np

def calculate_energies(prmtop_path, pdb_path, salt_conc_molar=0.15):
    print(f"Loading {prmtop_path} and {pdb_path}...")
    import parmed as pmd
    structure = pmd.load_file(prmtop_path)
    # We don't strictly need PDB positions for system creation, but we need them for context later.
    # Structure doesn't have positions unless loaded with inpcrd.
    # We will set positions on Context from the PDBFile.
    pdb = app.PDBFile(pdb_path)
    
    # Create system with standard settings (Microing core.py)
    print("Creating system via ParmEd...")
    system = structure.createSystem(nonbondedMethod=app.NoCutoff, 
                                 implicitSolvent=app.OBC2,
                                 implicitSolventSaltConc=salt_conc_molar*unit.molar)
    
    # Identify forces
    print("Forces in system:")
    nb_force = None
    custom_vdw_force = None
    gb_force = None
    
    for i, f in enumerate(system.getForces()):
        print(f"  {i}: {type(f).__name__}")
        if isinstance(f, openmm.NonbondedForce):
            nb_force = f
        elif isinstance(f, openmm.GBSAOBCForce):
            gb_force = f
        elif isinstance(f, openmm.CustomGBForce):
            gb_force = f
        elif isinstance(f, openmm.CustomNonbondedForce):
            # Assume first is VDW
            if custom_vdw_force is None:
                custom_vdw_force = f

    if custom_vdw_force:
        print(f"\nCustomNonbondedForce Exclusions: {custom_vdw_force.getNumExclusions()}")
        print(f"Energy Function: {custom_vdw_force.getEnergyFunction()}")
        
    if nb_force:
        print(f"NonbondedForce Exceptions: {nb_force.getNumExceptions()}")
        # Check EPSILON for first few particles to see if NonbondedForce includes VDW
        print("\nChecking NonbondedForce Particle Parameters (first 5):")
        for i in range(min(5, nb_force.getNumParticles())):
            chg, sig, eps = nb_force.getParticleParameters(i)
            print(f"  Atom {i}: Charge={chg}, Sigma={sig}, Epsilon={eps}")

        # Check Exceptions (1-4s)
        print("\nChecking NonbondedForce Exceptions (first 5):")
        for i in range(min(5, nb_force.getNumExceptions())):
            p1, p2, chg, sig, eps = nb_force.getExceptionParameters(i)
            print(f"  Exception {i}: Particles=({p1},{p2}), ChargeProd={chg}, Sigma={sig}, Epsilon={eps}")
    # Calculate initial VDW (Group 3)
    # Ensure CustomNonbondedForce is Group 3
    if custom_vdw_force:
        custom_vdw_force.setForceGroup(3)
        
    integrator = openmm.VerletIntegrator(0.001)
    context = openmm.Context(system, integrator)
    context.setPositions(pdb.positions)
    
    vdw_initial = context.getState(getEnergy=True, groups=(1<<3)).getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole)
    print(f"\nInitial VDW Energy (Group 3): {vdw_initial:.2f} kcal/mol")
    
    # Apply Fix Logic locally
    if nb_force and custom_vdw_force:
        print("\nApplying 1-4 VDW Correction (Test)...")
        
        vdw_14_force = openmm.CustomBondForce("4*epsilon*((sigma/r)^12 - (sigma/r)^6)")
        vdw_14_force.addPerBondParameter("sigma")
        vdw_14_force.addPerBondParameter("epsilon")
        vdw_14_force.setForceGroup(3) # Group 3
        
        count_moved = 0
        original_exceptions = []
        
        for i in range(nb_force.getNumExceptions()):
            p1, p2, chg, sig, eps = nb_force.getExceptionParameters(i)
            original_exceptions.append((p1, p2, chg, sig, eps))
            
            if eps._value != 0.0:
                # Add to CustomBondForce
                vdw_14_force.addBond(p1, p2, [sig, eps])
                # Exclude from CustomNonbondedForce (if not already)
                try:
                    custom_vdw_force.addExclusion(p1, p2)
                except Exception:
                    pass # Already excluded
                # Zero in NonbondedForce
                nb_force.setExceptionParameters(i, p1, p2, chg, sig, 0.0*unit.kilojoule_per_mole)
                count_moved += 1
                
        print(f"  Moved {count_moved} 1-4 interactions.")
        
        # Add force to system
        system.addForce(vdw_14_force)
        
        # Reinitialize context
        context.reinitialize(preserveState=True)
        
        vdw_corrected = context.getState(getEnergy=True, groups=(1<<3)).getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole)
        print(f"Corrected VDW Energy (Group 3): {vdw_corrected:.2f} kcal/mol")
        print(f"  (Should include CustomNonbondedForce + CustomBondForce)")
        
        # Breakdown
        # We need separate groups to see breakdown, but they are both in Group 3.
        # Let's temporarily move vdw_14 to Group 4
        vdw_14_force.setForceGroup(4)
        context.reinitialize(preserveState=True)
        vdw_main = context.getState(getEnergy=True, groups=(1<<3)).getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole)
        vdw_14 = context.getState(getEnergy=True, groups=(1<<4)).getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole)
        print(f"  Breakdown: Main VDW (excl 1-4) = {vdw_main:.2f}")
        print(f"  Breakdown: 1-4 VDW (Scaled)    = {vdw_14:.2f}")
        print(f"  Sum = {vdw_main + vdw_14:.2f}")
        
    return 0, 0, 0

if __name__ == "__main__":
    # Test on Complex Frame
    print("--- COMPLEX ---")
    calculate_energies('test/data/6t1h_1151_comp/complex.prmtop', 'test/results/6t1h_1151_comp/analysis_20260209_075721/temp_from_traj.pdb')

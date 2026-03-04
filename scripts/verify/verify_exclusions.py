
import openmm
import openmm.app as app
import sys

def verify_exclusions(prmtop_path):
    print(f"Loading {prmtop_path}...")
    prmtop = app.AmberPrmtopFile(prmtop_path)
    
    # Create system (Standard method, might not have CustomNonbondedForce unless we mimic core.py)
    # Be careful: standard createSystem puts everything in NonbondedForce.
    # We need to verify what happens in core.py's flow, BUT debugging the structure itself helps.
    # Actually, we need to inspect the System *as created by core.py*.
    # Using 'mmgbsa.core' to create system is better.
    
    # Let's fallback to standard Prmtop check first.
    # Creating system with standard options
    system_std = prmtop.createSystem(nonbondedMethod=app.NoCutoff)
    nb_std = [f for f in system_std.getForces() if isinstance(f, openmm.NonbondedForce)][0]
    
    # Count exceptions (1-2, 1-3, 1-4)
    n_exc_std = nb_std.getNumExceptions()
    print(f"Standard NonbondedForce Exceptions: {n_exc_std}")
    
    # Now simulate what happens if we had a CustomNonbondedForce
    # But wait, I can't check the LIVE system unless I load the one made by core.py? 
    # Or I can try to replicate the logic.
    
    # If I can't replicate easily, let's assume core.py uses ParmEd conversion or similar.
    # In core.py logs: "Replacing existing GB force (CustomGBForce)..."
    # This implies the input SYSTEM had CustomGBForce.
    
    # Let's try to load the system from the logic used in 'core.py' if possible, or just check standard behavior.
    # If standard behavior puts 1-4s in Exceptions, then they are Exclusions for CustomNonbondedForce IF copied.
    
    # Let's count 1-4s in Topology
    # 1-2, 1-3 are exclusions. 1-4 are scaled.
    # OpenMM Exceptions list usually contains ALL of them.
    # Exception: p1, p2, chgProd, sigma, epsilon.
    # If 1-2/1-3: params are 0,0,0.
    # If 1-4: params are scaled.
    
    count_14 = 0
    count_exclusion = 0
    for i in range(n_exc_std):
        p1, p2, c, s, e = nb_std.getExceptionParameters(i)
        if c._value == 0 and e._value == 0:
            count_exclusion += 1
        else:
            count_14 += 1
            
    print(f"  Exclusions (1-2/1-3): {count_exclusion}")
    print(f"  Exceptions (1-4): {count_14}")
    
    # If CustomNonbondedForce (VDW) is created, it needs `addExclusion(p1, p2)` for ALL of these?
    # Typical CustomNonbondedForce logic:
    # "interaction is computed for every pair... excluding those in the recommended exclusion list"
    # The exclusion list SHOULD include 1-2, 1-3 AND 1-4 if we want to handle 1-4s separately or scaled.
    # If we DON'T exclude 1-4s, CustomNonbondedForce computes Full VDW for them.
    
    print("\nHypothesis Check:")
    print(f"If CustomNonbondedForce only excludes 1-2/1-3 ({count_exclusion}),")
    print(f"Then {count_14} 1-4 pairs are computed as FULL VDW interactions.")
    
    return count_14

if __name__ == "__main__":
    verify_exclusions('test/data/6t1h_1151_comp/complex.prmtop')

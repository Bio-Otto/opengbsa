
import parmed as pmd
from openmm import app, openmm, unit
import sys
import time

def verify_workaround():
    print("Starting verify_workaround...", flush=True)
    prmtop = "test/data/6t1h_6466_comp/complex.prmtop"
    
    print(f"Loading {prmtop}...", flush=True)
    struct = pmd.load_file(prmtop)
    
    # Simulate the logic in core.py
    nonbondedCutoff = 12.0*unit.angstroms
    nonbondedMethod = app.CutoffNonPeriodic
    implicitSolvent = app.OBC2
    implicitSolventSaltConc = 0.15*unit.molar
    
    print("\n--- Applying Workaround Logic ---", flush=True)
    
    # WORKAROUND START
    pass_implicit = implicitSolvent
    pass_salt = implicitSolventSaltConc
    
    if nonbondedCutoff is not None and implicitSolvent is not None:
         print("Detected Cutoff + ImplicitSolvent. Disabling ImplicitSolvent for createSystem to avoid hang.", flush=True)
         pass_implicit = None
         pass_salt = None
    # WORKAROUND END
    
    print(f"Calling createSystem with implicitSolvent={pass_implicit}...", flush=True)
    start = time.time()
    try:
        sys = struct.createSystem(nonbondedMethod=nonbondedMethod,
                                  nonbondedCutoff=nonbondedCutoff,
                                  implicitSolvent=pass_implicit,
                                  implicitSolventSaltConc=pass_salt)
        print(f"Success! Time: {time.time()-start:.2f} s", flush=True)
    except Exception as e:
        print(f"Failed: {e}", flush=True)
        return

    # Phase 2: Refine Forces (Manual Addition)
    if pass_implicit is None and implicitSolvent is not None:
        print("Now simulating GBSA refinement...", flush=True)
        # Mocking refinement
        gb_force = openmm.GBSAOBCForce()
        gb_force.setNonbondedMethod(openmm.GBSAOBCForce.NoCutoff) # Implicit solvent is always NoCutoff
        # ... (rest of simple addition)
        sys.addForce(gb_force)
        print("Added GBSA force manually.", flush=True)

    print("\nVerification Complete.", flush=True)

if __name__ == "__main__":
    verify_workaround()

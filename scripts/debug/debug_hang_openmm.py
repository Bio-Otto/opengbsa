
from openmm import app, openmm, unit
import sys
import time

def debug_hang_openmm():
    print("Starting debug_hang_openmm...")
    prmtop_path = "test/data/6t1h_6466_comp/complex.prmtop"
    
    print(f"Loading {prmtop_path}...")
    prmtop = app.AmberPrmtopFile(prmtop_path)
    
    # Test 1: Cutoff + Implicit Solvent
    print("\n--- Test 1: AmberPrmtopFile.createSystem (Cutoff=12.0 + OBC2) ---")
    start = time.time()
    try:
        sys = prmtop.createSystem(nonbondedMethod=app.CutoffNonPeriodic, 
                                  nonbondedCutoff=12.0*unit.angstroms,
                                  implicitSolvent=app.OBC2)
        print(f"Success! Time: {time.time()-start:.2f} s")
    except Exception as e:
        print(f"Failed: {e}")

    print("\nDebug Complete.")

if __name__ == "__main__":
    debug_hang_openmm()


import parmed as pmd
from openmm import app, openmm, unit
import sys
import time

def debug_hang():
    print("Starting debug_hang...")
    prmtop = "test/data/6t1h_6466_comp/complex.prmtop"
    
    print(f"Loading {prmtop}...")
    struct = pmd.load_file(prmtop)
    print(f"Loaded {len(struct.atoms)} atoms.")
    
    # Test 1: NoCutoff (Standard)
    print("\n--- Test 1: createSystem (NoCutoff) ---")
    start = time.time()
    try:
        sys = struct.createSystem(nonbondedMethod=app.NoCutoff)
        print(f"Success! Time: {time.time()-start:.2f} s")
    except Exception as e:
        print(f"Failed: {e}")

    # Test 3: Cutoff + Implicit Solvent
    print("\n--- Test 3: createSystem (Cutoff=12.0 + OBC2) ---")
    start = time.time()
    try:
        sys = struct.createSystem(nonbondedMethod=app.CutoffNonPeriodic, 
                                  nonbondedCutoff=12.0*unit.angstroms,
                                  implicitSolvent=app.OBC2)
        print(f"Success! Time: {time.time()-start:.2f} s")
    except Exception as e:
        print(f"Failed: {e}")

    print("\nDebug Complete.")

if __name__ == "__main__":
    debug_hang()

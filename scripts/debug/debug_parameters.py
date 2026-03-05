
import parmed as pmd
import sys

def check_parameters(prmtop_path):
    print(f"Loading {prmtop_path}...")
    struct = pmd.load_file(prmtop_path)
    
    print(f"Structure type: {type(struct)}")
    
    # Check first few atoms
    print("\nChecking atom parameters (epsilon, sigma):")
    zero_eps = 0
    total = 0
    for atom in struct.atoms:
        total += 1
        if atom.epsilon == 0:
            zero_eps += 1
            
    print(f"Total atoms: {total}")
    print(f"Atoms with epsilon=0: {zero_eps}")
    
    if zero_eps > 0:
        print("  -> CONFIRMED: Atoms have zero epsilon.")
        
    # Check if LJ_depth exists
    if hasattr(struct, 'LJ_depth'):
        print("\nAmberParm has LJ_depth array. Checking values...")
        non_zero_depth = sum(1 for x in struct.LJ_depth if x > 0)
        print(f"Entries in LJ_depth > 0: {non_zero_depth}")
        if non_zero_depth > 0 and zero_eps == total:
             print("  -> DISCREPANCY: Atom.epsilon is 0, but LJ_depth has values. Need to sync!")

if __name__ == "__main__":
    check_parameters("test/data/6t1h_6466_comp/complex.prmtop")

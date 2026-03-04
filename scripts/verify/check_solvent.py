
import parmed as pmd
import sys

def check_solvent(prmtop_path):
    print(f"Checking {prmtop_path}...")
    struct = pmd.load_file(prmtop_path)
    print(f"Total Atoms: {len(struct.atoms)}")
    print(f"Total Residues: {len(struct.residues)}")
    
    solvent_residues = ['WAT', 'HOH', 'TIP3', 'SOL', 'Na+', 'Cl-', 'NA', 'CL', 'K+', 'Mg+']
    
    found_solvent = []
    for res in struct.residues:
        if res.name in solvent_residues:
            found_solvent.append(res.name)
            
    if found_solvent:
        from collections import Counter
        counts = Counter(found_solvent)
        print("CRITICAL: Found solvent/ion residues:")
        for name, count in counts.items():
            print(f"  {name}: {count}")
    else:
        print("No standard solvent/ion residues found.")

if __name__ == "__main__":
    check_solvent("test/data/6t1h_6466_comp/complex.prmtop")

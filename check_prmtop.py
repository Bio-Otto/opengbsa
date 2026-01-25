
import parmed as pmd
import sys

try:
    path = '/home/bio-otto/Desktop/OXA-MD/5iy2_6466/replica_1/topol_dry.prmtop'
    print(f"Loading {path}...")
    mol = pmd.load_file(path)
    
    print("Checking first 5 atoms:")
    for i in range(5):
        atom = mol.atoms[i]
        print(f"Atom {i} ({atom.name}): sigma={atom.sigma}, epsilon={atom.epsilon}, rmin={atom.rmin}, type={atom.type}")
        
    print("\nChecking Ligand atoms:")
    lig = [a for a in mol.atoms if a.residue.name == 'LIG']
    if lig:
        for i in range(min(5, len(lig))):
             atom = lig[i]
             print(f"Lig Atom {i} ({atom.name}): sigma={atom.sigma}, epsilon={atom.epsilon}")
    else:
        print("No LIG found")

except Exception as e:
    print(e)

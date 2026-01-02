
import mdtraj as md
import sys

def check_termini(gro_file):
    t = md.load(gro_file)
    prot = t.topology.subset(t.topology.select("protein"))
    
    res = list(prot.residues)
    first = res[0]
    last = res[-1]
    
    print(f"First Residue: {first} (Index {first.index})")
    print(f"Last Residue: {last} (Index {last.index})")
    
    # Check if PRO265 and THR25 match these
    # Note: Indices might shift in subsetting
    
    print("Checking Bonds involving termini...")
    for bond in prot.bonds:
        r1 = bond[0].residue
        r2 = bond[1].residue
        
        # Check for long range bonds (index diff > 1)
        if abs(r1.index - r2.index) > 1:
            print(f"Found long-range bond: {r1}-{r1.index} -- {r2}-{r2.index}")
            print(f"  Atoms: {bond[0]} -- {bond[1]}")

if __name__ == "__main__":
    check_termini(sys.argv[1])


import mdtraj as md
import sys

def check_bonds(gro_file):
    print(f"Loading GRO: {gro_file}")
    t = md.load(gro_file)
    
    # Select dry
    dry = t.topology.select("protein or resname UNL or resname LIG")
    t_dry = t.atom_slice(dry)
    
    n_bonds = 0
    # MDTraj topology interaction
    # bonds is an iterator
    bonds = list(t_dry.topology.bonds)
    n_bonds = len(bonds)
    
    print(f"Number of atoms: {t_dry.n_atoms}")
    print(f"Number of bonds: {n_bonds}")
    
    if n_bonds == 0:
        print("❌ NO BONDS FOUND! This is why image_molecules failed.")
        print("MDTraj standard bond connect did not work on this GRO subset.")
    else:
        print(f"✅ Bonds found: {n_bonds}. MDTraj connectivity seems OK.")

if __name__ == "__main__":
    check_bonds(sys.argv[1])

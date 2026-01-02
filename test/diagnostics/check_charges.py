
import sys
import mdtraj as md
import numpy as np
from openmm.app import PDBFile

def check_structure_charge(pdb_file):
    print(f"Checking charges for: {pdb_file}")
    
    pdb = PDBFile(pdb_file)
    topology = pdb.topology
    
    total_charge = 0.0
    ligand_charge = 0.0
    protein_charge = 0.0
    ligand_resname = 'UNL' # Assumption
    
    # Try to find ligand resname
    resnames = set([r.name for r in topology.residues()])
    print(f"Residues found: {len(resnames)} types")
    if 'UNL' not in resnames:
        print(f"Warning: UNL not found. Available: {list(resnames)[:10]}...")
        # Simple heuristic
        for r in resnames:
            if len(r) == 3 and r not in ['ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL', 'HOH', 'WAT', 'NA', 'CL']:
                ligand_resname = r
                break
    print(f"Assuming Ligand Resname: {ligand_resname}")

    # PDB files don't usually have full partial charges unless prepared by Amber/Charmm
    # But we can look at formal charges if present, or guess based on atom names/residues
    
    print("\n--- Atom Charge Analysis (from PDB columns if available) ---")
    # OpenMM PDBFile does not readily expose partial charges unless in specific columns
    # Let's use simple residue counting for protein
    
    pos_res = ['ARG', 'LYS', 'HIP']
    neg_res = ['ASP', 'GLU']
    
    net_prot = 0
    for res in topology.residues():
        if res.name in pos_res: net_prot += 1
        elif res.name in neg_res: net_prot -= 1
        
    print(f"Estimated Protein Net Charge (at pH 7): {net_prot}")
    print("  (Based on standard ARG/LYS +1 and ASP/GLU -1 counting)")
    
    print("\n--- Ligand AM1-BCC Fallback Check ---")
    print("The previous error indicated AM1-BCC failing. This usually happens if:")
    print("1. Molecules have strange valency.")
    print("2. Required toolkits (AmberTools/OpenEye) are missing.")
    print("3. The ligand structure (SDF) has disconnected fragments.")
    
if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python check_charges.py complex.pdb")
    else:
        check_structure_charge(sys.argv[1])

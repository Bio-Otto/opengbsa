
import mdtraj as md

input_pdb = 'test/3if6_test/solv_ions.pdb'
output_pdb = 'test/3if6_test/solv_ions_split.pdb'

print(f"Loading topology from {input_pdb}...")
traj = md.load(input_pdb)

# Get the PDB residue number (resSeq) of index 242
# Note: MDTraj 0-indexed residues map to whatever is in the PDB
res_242 = traj.topology.residue(242)
res_Seq_start = res_242.resSeq

print(f"Residue Index 242 corresponds to PDB Residue Number: {res_Seq_start} (Name: {res_242.name})")

# Read PDB and rewrite chains
print("Rewriting PDB chains...")
with open(input_pdb, 'r') as f:
    lines = f.readlines()

with open(output_pdb, 'w') as f:
    current_chain = 'A'
    for line in lines:
        if line.startswith("ATOM") or line.startswith("HETATM"):
            try:
                res_seq = int(line[22:26])
                # If we hit the residue number corresponding to the split point
                # Make sure we don't switch back if numbers wrap, but usually they increase.
                # Assuming sequential.
                
                # Checking PDB residue number directly (from grep inspection)
                if res_seq >= 242:
                    current_chain = 'B'
                else:
                    current_chain = 'A'
                
                # Write line with new chain ID (Col 21 in 0-index -> Col 22 in 1-index)
                # PDB Format: Chain ID is column 21 (0-based)
                new_line = line[:21] + current_chain + line[22:]
                f.write(new_line)
            except:
                f.write(line)
        else:
            f.write(line)

print(f"Saved split PDB to {output_pdb}")

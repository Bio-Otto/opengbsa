import sys

input_file = "test/data/3if6_dimer_test/3if6.pdb"
output_file = "test/data/3if6_dimer_test/3if6_fixed.pdb"

lines = []
with open(input_file) as f:
    for line in f:
        if line.startswith("ATOM  ") or line.startswith("HETATM"):
            resname = line[17:20].strip()
            if resname in ["SOL", "WAT", "HOH", "NA", "CL", "K", "MG", "CA", "UNL"]:
                # skip solvent and ligand; we just want the protein to be properly sequential
                continue
            lines.append(line)

chains = [] # list of lists of lines
current_res_atoms = [] # sets of atom names seen for each chain

last_resid = None

for line in lines:
    resid = line[22:27].strip() # include insertion code if any
    atom_name = line[12:16].strip()
    alt_loc = line[16]
    
    # We only care about alt_loc A or ' '
    if alt_loc not in [' ', 'A']:
        continue
        
    identifier = atom_name
    
    if resid != last_resid:
        # new residue block starts, reset seen atom sets
        current_res_atoms = [set() for _ in range(max(1, len(chains)))]
        last_resid = resid
         
    # find which chain to put it in
    placed = False
    for i, seen in enumerate(current_res_atoms):
        if identifier not in seen:
            seen.add(identifier)
            while len(chains) <= i:
                chains.append([])
            chains[i].append(line)
            placed = True
            break
            
    if not placed:
        # Need a new chain instance for this duplicated residue
        current_res_atoms.append({identifier})
        chains.append([line])

print(f"Discovered {len(chains)} interwoven chains!")

# Write them out
with open(output_file, "w") as f:
    for i, chain_lines in enumerate(chains):
        chain_id = chr(ord('A') + i) if i < 26 else 'X'
        for line in chain_lines:
            new_line = line[:21] + chain_id + line[22:]
            f.write(new_line)
        f.write("TER\n")
        print(f"Chain {chain_id} wrote {len(chain_lines)} atoms.")
        
print("Saved to 3if6_fixed.pdb")

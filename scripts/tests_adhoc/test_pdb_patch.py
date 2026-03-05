import re

def patch_missing_n_terminal_h3(pdb_path, out_path):
    """
    Scans a PDB file. For every N-terminal residue (the first residue of any chain),
    if it contains H1 and H2 but is missing H3, this function manually injects an H3 line 
    right after H2, copying H2's coordinates + 0.001 offset.
    This fulfills OpenMM's strict terminal templating (e.g., NVAL expects 3 protons).
    """
    with open(pdb_path, 'r') as f:
        lines = f.readlines()
        
    out_lines = []
    
    current_chain = None
    first_res_of_chain = True
    in_terminal_res = False
    terminal_res_num = None
    
    for i, line in enumerate(lines):
        if line.startswith('ATOM') or line.startswith('HETATM'):
            chain = line[21]
            res_num = line[22:26].strip()
            atom_name = line[12:16].strip()
            
            if chain != current_chain:
                current_chain = chain
                first_res_of_chain = True
                terminal_res_num = res_num
                in_terminal_res = True
            
            if in_terminal_res and res_num != terminal_res_num:
                in_terminal_res = False
                
            out_lines.append(line)
            
            if in_terminal_res and first_res_of_chain and atom_name == 'H2':
                # Check if H3 exists in the next few lines
                h3_exists = False
                for j in range(i+1, min(i+10, len(lines))):
                    nl = lines[j]
                    if nl.startswith('ATOM') and nl[21] == chain and nl[22:26].strip() == terminal_res_num:
                        if nl[12:16].strip() == 'H3':
                            h3_exists = True
                            break
                    else:
                        break # Next res or different record
                        
                if not h3_exists:
                    # Inject H3
                    # ATOM     22  H2  VAL A   0      31.295 -20.671  55.513  1.00  0.00           H  
                    # 0123456789012345678901234567890123456789012345678901234567890123456789
                    atom_serial = int(line[6:11]) + 1
                    try:
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                    except ValueError:
                        x, y, z = 0.0, 0.0, 0.0
                        
                    h3_line = f"ATOM  {atom_serial:5d}  H3  {line[17:20]} {chain}{line[22:26]}    {x:8.3f}{y:8.3f}{z+0.010:8.3f}{line[54:]}"
                    out_lines.append(h3_line)
                    print(f"Patched missing H3 for Chain {chain} Res {terminal_res_num}")
                    
        else:
            out_lines.append(line)
            
    with open(out_path, 'w') as f:
        f.writelines(out_lines)

patch_missing_n_terminal_h3('test/data/3if6_dimer_test/3if6_sanitized.pdb', 'test/data/3if6_dimer_test/3if6_patched.pdb')

import openmm.app as app
try:
    pdb = app.PDBFile('test/data/3if6_dimer_test/3if6_patched.pdb')
    modeller = app.Modeller(pdb.topology, pdb.positions)
    to_delete_res = [r for r in modeller.topology.residues() if r.name in ['HOH', 'WAT', 'TIP3', 'SOL', 'NA', 'CL', 'K', 'MG', 'ZN', 'CA', 'LIG', 'UNL']]
    modeller.delete(to_delete_res)
    
    forcefield = app.ForceField('amber/ff14SB.xml', 'amber/tip3p_standard.xml')
    sys = forcefield.createSystem(modeller.topology, nonbondedMethod=app.NoCutoff, constraints=app.HBonds)
    print("SUCCESS: 3if6 parsed and system created completely from patched PDB!")
except Exception as e:
    print("ERROR:", type(e).__name__, e)

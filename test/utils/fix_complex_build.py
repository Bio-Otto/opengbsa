#!/usr/bin/env python3
"""Fix build_complex_system method"""

# Read the file
with open('mmgbsa/core.py', 'r') as f:
    lines = f.readlines()

# Find the section to replace (around line 1230-1242)
# We need to replace the topology building part in build_complex_system

new_section = '''        # Extract only protein atoms (exclude ligand) and group by residue
        protein_atoms = [atom for atom in protein_pdbfile.topology.atoms() if atom.residue.name != 'LIG']
        
        # Group atoms by residue to ensure contiguous residues
        from collections import defaultdict
        residue_atoms = defaultdict(list)
        for atom in protein_atoms:
            residue_atoms[(atom.residue.chain.id, atom.residue.id)].append(atom)
        
        # Sort residues and atoms within each residue
        sorted_residue_keys = sorted(residue_atoms.keys())
        protein_atoms_sorted = []
        for res_key in sorted_residue_keys:
            atoms_in_residue = sorted(residue_atoms[res_key], key=lambda a: a.index)
            protein_atoms_sorted.extend(atoms_in_residue)
        
        protein_top = app.Topology()
        chain_map = {}
        residue_map = {}
        atom_map = {}
        
        for atom in protein_atoms_sorted:
            chain = atom.residue.chain
            if chain.id not in chain_map:
                chain_map[chain.id] = protein_top.addChain(chain.id)
            
            res_key = (chain.id, atom.residue.id)
            if res_key not in residue_map:
                residue_map[res_key] = protein_top.addResidue(atom.residue.name, chain_map[chain.id], id=atom.residue.id)
            
            new_atom = protein_top.addAtom(atom.name, atom.element, residue_map[res_key], id=atom.id)
            atom_map[atom] = new_atom
        
        # Add bonds for protein using atom_map
        for bond in protein_pdbfile.topology.bonds():
            if bond[0].residue.name != 'LIG' and bond[1].residue.name != 'LIG':
                if bond[0] in atom_map and bond[1] in atom_map:
                    try:
                        protein_top.addBond(atom_map[bond[0]], atom_map[bond[1]])
                    except Exception:
                        continue
        
        # Extract protein positions in sorted order
        all_positions = protein_pdbfile.positions.value_in_unit(unit.nanometer)
        protein_coords = [all_positions[atom.index] for atom in protein_atoms_sorted]
        protein_positions = unit.Quantity(protein_coords, unit.nanometer)
'''

# Find and replace lines 1230-1260
start_line = 1229  # 0-indexed, so line 1230
end_line = 1260

# Replace the section
new_lines = lines[:start_line] + [new_section] + lines[end_line:]

# Write back
with open('mmgbsa/core.py', 'w') as f:
    f.writelines(new_lines)

print(f'✓ Replaced lines {start_line+1}-{end_line} in core.py')

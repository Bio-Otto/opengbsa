import openmm.app as app
import openmm.unit as unit

pdb = app.PDBFile('test/data/3if6_dimer_test/3if6_sanitized.pdb')

modeller = app.Modeller(pdb.topology, pdb.positions)
to_delete_res = [r for r in modeller.topology.residues() if r.name in ['HOH', 'WAT', 'TIP3', 'SOL', 'NA', 'CL', 'K', 'MG', 'ZN', 'CA', 'LIG', 'UNL']]
modeller.delete(to_delete_res)

# Let's manually add the missing H3 atom to Residue 0 (VAL)
val_res = next(r for r in modeller.topology.residues() if r.index == 0)

# Check existing hydrogens in VAL 0
existing_h = [a.name for a in val_res.atoms() if a.element.symbol == 'H']
print("Existing H in VAL 0:", existing_h)

# The Amber NVAL template expects H1, H2, H3 attached to N.
# The PDB only has H1, H2. We need to add H3.
n_atom = next(a for a in val_res.atoms() if a.name == 'N')

# Add H3 atom to topology
h_element = app.Element.getBySymbol('H')
h3_atom = modeller.topology.addAtom('H3', h_element, val_res)
modeller.topology.addBond(n_atom, h3_atom)

# Give H3 a position (just copy N position and offset slightly)
import numpy as np
n_idx = n_atom.index
n_pos = modeller.positions[n_idx].value_in_unit(unit.nanometers)

# Append position
# Position is a list of Vec3, we need to append one more
positions_list = list(modeller.positions.value_in_unit(unit.nanometers))
positions_list.append(n_pos + np.array([0.1, 0.1, 0.1])) 
modeller.positions = positions_list * unit.nanometers

# Check ASN 241
asn_res = next(r for r in modeller.topology.residues() if r.index == 241)
# Does it have OXT?
has_oxt = any(a.name == 'OXT' for a in asn_res.atoms())
# In 'test_h2.py', ASN failed `addHydrogens`. Will it fail `createSystem` without it?
print("ASN 241 has OXT?", has_oxt)

forcefield = app.ForceField('amber/ff14SB.xml', 'amber/tip3p_standard.xml')
try:
    sys = forcefield.createSystem(modeller.topology, nonbondedMethod=app.NoCutoff, constraints=app.HBonds)
    print("System created successfully!")
except Exception as e:
    print("create_system ERROR:", type(e).__name__, e)

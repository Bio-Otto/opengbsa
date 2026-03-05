import openmm.app as app

pdb = app.PDBFile('test/data/3if6_dimer_test/3if6_patched.pdb')
modeller = app.Modeller(pdb.topology, pdb.positions)

# Let's see what ASN 241 is bonded to
casn = None
for r in modeller.topology.residues():
    if r.index == 241:
        casn = r
        break

print(f"Bonds for {casn.name} {casn.index}:")
for bond in modeller.topology.bonds():
    a1, a2 = bond
    if (a1.residue == casn and a2.residue != casn) or (a2.residue == casn and a1.residue != casn):
        other = a2 if a1.residue == casn else a1
        my_atom = a1 if a1.residue == casn else a2
        print(f"  {my_atom.name} is bonded to EXT: {other.residue.name} {other.residue.index} atom {other.name}")


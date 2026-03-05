import parmed as pmd

struct = pmd.load_file('test/data/3if6_dimer_test/3if6.pdb')

print("Generating clean PDB using ParmEd...")

# ParmEd usually assigns distinct residue numbers if it detects breaks in bonds or chain splits
# Let's save it directly to see if ParmEd fixes the interwoven residues
struct.save('test/data/3if6_dimer_test/3if6_cleaned.pdb', overwrite=True)
print("Saved 3if6_cleaned.pdb")

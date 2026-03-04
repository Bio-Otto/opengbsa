import MDAnalysis as mda
import warnings

warnings.filterwarnings('ignore')

u = mda.Universe('test/data/3if6_dimer_test/3if6.pdb')

print("Generating segments using MDAnalysis...")
# MDAnalysis guesses bonds and can split based on fragments automatically.
segments = u.atoms.fragments
print(f"Found {len(segments)} fragments.")

# Assign unique segment IDs to each fragment (which become chain IDs in PDB)
letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
for i, frag in enumerate(segments):
    if i < len(letters):
        frag.segments.segids = letters[i]
    else:
        frag.segments.segids = f"X{i}"

# Also, ensure unique residue numbers for waters
waters = u.select_atoms('resname SOL or resname HOH or resname WAT or resname TIP3')
if len(waters) > 0:
    print(f"Found {len(waters.residues)} unique water residues.")
    # Assign unique resids sequentially
    for i, res in enumerate(waters.residues):
        res.resid = 9999 + i + 1

# Save it
u.atoms.write('test/data/3if6_dimer_test/3if6_cleaned_mda.pdb')
print("Saved 3if6_cleaned_mda.pdb")

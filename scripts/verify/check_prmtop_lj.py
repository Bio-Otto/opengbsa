import parmed as pmd
print("Loading complex.prmtop...")
struct = pmd.load_file('test/data/6t1h_1151_comp/complex.prmtop')
print(f"Loaded {len(struct.atoms)} atoms.")

non_zero_eps = 0
for i, atom in enumerate(struct.atoms):
    if atom.epsilon != 0.0:
        non_zero_eps += 1
    if i < 10:
        print(f"Atom {i} ({atom.name}): Eps={atom.epsilon}, Rmin={atom.rmin}")

print(f"Total atoms with Non-Zero Epsilon: {non_zero_eps}")

import mdtraj as md
import os

prmtop = 'test/data/6t1h_prmtop_test/complex.prmtop'
gro = 'test/data/6t1h_prmtop_test/md_complex_prod.gro'
xtc = 'test/data/6t1h_prmtop_test/md_complex_prod.xtc'

print(f"Checking {prmtop}...")
try:
    t1 = md.load_prmtop(prmtop)
    print(f"Prmtop Atoms: {t1.n_atoms}")
except Exception as e:
    print(f"Prmtop Error: {e}")

print(f"Checking {gro}...")
try:
    t2 = md.load(gro)
    print(f"GRO Atoms: {t2.n_atoms}")
except Exception as e:
    print(f"GRO Error: {e}")

print(f"Checking XTC with GRO...")
try:
    t3 = md.load(xtc, top=gro)
    print(f"XTC (with GRO) Atoms: {t3.n_atoms}, Frames: {t3.n_frames}")
except Exception as e:
    print(f"XTC Error: {e}")

print(f"Checking XTC with PRMTOP...")
try:
    t4 = md.load(xtc, top=prmtop)
    print(f"XTC (with PRMTOP) Atoms: {t4.n_atoms}, Frames: {t4.n_frames}")
except Exception as e:
    print(f"XTC Error: {e}")

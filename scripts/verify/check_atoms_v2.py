import mdtraj as md
import traceback

prmtop = 'test/data/6t1h_prmtop_test/complex.prmtop'
gro = 'test/data/6t1h_prmtop_test/md_complex_prod.gro'
xtc = 'test/data/6t1h_prmtop_test/md_complex_prod.xtc'

print("--- DEBUG ---")
try:
    md.load(xtc, top=gro)
    print("Success loading with GRO")
except Exception:
    traceback.print_exc()

print("\n--- DEBUG PRMTOP ---")
try:
    md.load(xtc, top=prmtop)
    print("Success loading with PRMTOP")
except Exception:
    traceback.print_exc()

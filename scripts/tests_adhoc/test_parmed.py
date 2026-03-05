
import sys
print("Starting test_parmed...", flush=True)

try:
    import parmed as pmd
    print("Imported parmed successfully.", flush=True)
except Exception as e:
    print(f"Failed to import parmed: {e}", flush=True)

print("Attempting to load file...", flush=True)
try:
    s = pmd.load_file("test/data/6t1h_6466_comp/complex.prmtop")
    print(f"Loaded file with {len(s.atoms)} atoms.", flush=True)
except Exception as e:
    print(f"Failed to load file: {e}", flush=True)

print("Test complete.", flush=True)

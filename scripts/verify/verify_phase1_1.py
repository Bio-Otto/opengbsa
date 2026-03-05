
import pandas as pd
import glob
import os

# Find most recent results
files = glob.glob('test/results/6t1h_1151_comp/analysis_*/fixed_enhanced_mmgbsa_results_gbn.csv')
if not files:
    print("No results found!")
    exit(1)
latest = max(files, key=os.path.getctime)

print(f"Analyzing most recent result: {latest}")
df = pd.read_csv(latest)

print("="*60)
print("PHASE 1.1 VERIFICATION (After VDW Fix)")
print("="*60)

binding = df['binding_energy'].mean()
vdw = df['delta_vdw'].mean()
ele = df['delta_elec'].mean()
gb = df['delta_gb'].mean()
sa = df['delta_sa'].mean()

print(f"\nBinding Energy: {binding:.2f} ± {df['binding_energy'].std():.2f}")
print(f"\nComponent Breakdown:")
print(f"  VDW:  {vdw:.2f} ± {df['delta_vdw'].std():.2f}")
print(f"  ELE:  {ele:.2f} ± {df['delta_elec'].std():.2f}")
print(f"  GB:   {gb:.2f} ± {df['delta_gb'].std():.2f}")
print(f"  SA:   {sa:.2f} ± {df['delta_sa'].std():.2f}")

print(f"\n{'='*60}")
print("COMPARISON WITH MMPBSA.py (Reference)")
print("="*60)
print(f"{'Component':<15} {'OpenGBSA':<12} {'MMPBSA.py':<12} {'Diff':<12}")
print("-"*60)
# MMPBSA.py Reference Values
ref_binding = -11.67
ref_vdw = -19.27
ref_ele = -6.39
ref_gb = 16.67
ref_sa = -2.68

print(f"{'Binding':<15} {binding:<12.2f} {ref_binding:<12.2f} {binding - ref_binding:<12.2f}")
print(f"{'VDW':<15} {vdw:<12.2f} {ref_vdw:<12.2f} {vdw - ref_vdw:<12.2f}")
print(f"{'ELE':<15} {ele:<12.2f} {ref_ele:<12.2f} {ele - ref_ele:<12.2f}")
print(f"{'GB':<15} {gb:<12.2f} {ref_gb:<12.2f} {gb - ref_gb:<12.2f}")
print(f"{'SA':<15} {sa:<12.2f} {ref_sa:<12.2f} {sa - ref_sa:<12.2f}")

import pandas as pd
import os
import glob

# Find latest analysis
results_dir = "test/results/6t1h_prmtop_analysis"
list_of_dirs = glob.glob(os.path.join(results_dir, "analysis_*"))
latest_dir = max(list_of_dirs, key=os.path.getctime)

print(f"Checking results in: {latest_dir}")

complex_csv = os.path.join(latest_dir, "per_residue_complex.csv")
receptor_csv = os.path.join(latest_dir, "per_residue_receptor.csv")
ligand_csv = os.path.join(latest_dir, "per_residue_ligand.csv")

if not os.path.exists(complex_csv):
    print("❌ per_residue_complex.csv NOT found.")
    exit(1)

print("✅ per_residue_complex.csv found.")

df = pd.read_csv(complex_csv)
print(f"Columns: {df.columns.tolist()}")

if 'complex_total' not in df.columns:
    print("❌ complex_total column missing.")
    exit(1)

mean_complex = df['complex_total'].mean()
print(f"Mean Complex Total Energy: {mean_complex:.2f} kcal/mol")
# This should be large negative (Interaction with protein + solvation).
# Note: "Complex Total" includes Solvation. Large negative.
# Interaction (VdW + Ele) is also large negative (attractive).

top = df.sort_values('complex_total').head(5)
print("\nTop 5 Most Favorable Residues (Complex Total):")
print(top[['residue_name','residue_number','complex_total', 'complex_vdw', 'complex_electrostatic', 'complex_solvation']])

# Check Ligand CSV (Standard Binding)
if os.path.exists(ligand_csv):
    df_lig = pd.read_csv(ligand_csv)
    print("\nTop 5 Most Favorable Residues (Standard Binding):")
    print(df_lig.sort_values('total').head(5)[['residue_name','residue_number','total']])

print("\nVerification Complete.")


import pandas as pd
import numpy as np

def load_mmpbsa_dat(filepath):
    """
    Parses MMPBSA.py decomposition file manually since it's a fixed-width/CSV hybrid.
    Extracts residue-level total and component energies.
    """
    data = []
    with open(filepath, 'r') as f:
        lines = f.readlines()
    
    # Skip header lines until "Residue,Internal"
    start_idx = 0
    for i, line in enumerate(lines):
        if line.strip().startswith("Residue,Internal"):
            start_idx = i + 2 # Skip header and unit line
            break
            
    for line in lines[start_idx:]:
        parts = line.strip().split(',')
        if len(parts) < 2: continue
        
        # Parse Residue Name and Number (e.g. "SER   1")
        res_str = parts[0]
        res_name = res_str[:3]
        try:
            res_num = int(res_str[3:].strip())
        except ValueError:
            continue
            
        # MMPBSA.py output structure based on header in file:
        # Residue, Internal, VDW, Electrostatic, Polar Solvation, Non-Polar Solv, TOTAL
        # Each has Avg, Std Dev, Std Err.
        # Indices (0-based from split):
        # 0: Residue
        # 1-3: Internal (Avg, SD, SE)
        # 4-6: VDW
        # 7-9: Electrostatic
        # 10-12: Polar Solvation
        # 13-15: Non-Polar Solvation
        # 16-18: TOTAL
        
        try:
            vdw = float(parts[4])
            ele = float(parts[7])
            gb = float(parts[10]) # Polar
            sa = float(parts[13]) # Non-Polar
            total = float(parts[16])
            
            data.append({
                'residue_number': res_num,
                'residue_name': res_name,
                'mmpbsa_vdw': vdw,
                'mmpbsa_ele': ele,
                'mmpbsa_gb': gb,
                'mmpbsa_sa': sa,
                'mmpbsa_total': total
            })
        except (ValueError, IndexError):
            continue
            
    return pd.DataFrame(data)

# 1. Global Analysis
print("--- GLOBAL ENERGY COMPARISON (kcal/mol) ---")
try:
    opengbsa_df = pd.read_csv('test/results/6t1h_1151_comp/analysis_20260207_223534/fixed_enhanced_mmgbsa_results_gbn.csv')
    
    # Calculate means
    components = {
        'VDW': opengbsa_df['delta_vdw'].mean(),
        'EEL': opengbsa_df['delta_elec'].mean(),
        'GB': opengbsa_df['delta_gb'].mean(),
        'SA': opengbsa_df['delta_sa'].mean(),
        'Total': opengbsa_df['binding_energy'].mean()
    }
    
    print(f"{'Component':<10} {'OpenGBSA':<15} {'MMPBSA.py (Ref)':<15} {'Diff':<15}")
    mmpbsa_global = {
        'VDW': -19.2711,
        'EEL': -6.3924,
        'GB': 16.6683,
        'SA': -2.6771,
        'Total': -11.6724
    }
    
    for comp in ['VDW', 'EEL', 'GB', 'SA', 'Total']:
        og_val = components.get(comp, 0.0)
        ref_val = mmpbsa_global.get(comp, 0.0)
        print(f"{comp:<10} {og_val:<15.4f} {ref_val:<15.4f} {og_val - ref_val:<15.4f}")

except Exception as e:
    print(f"Error analyzing global results: {e}")

# 2. Per-Residue Analysis
print("\n--- PER-RESIDUE COMPARISON (Top Differences) ---")
try:
    # Load OpenGBSA
    og_res_df = pd.read_csv('test/results/6t1h_1151_comp/analysis_20260207_223534/per_residue_detailed.csv')
    
    # Load MMPBSA
    mmpbsa_res_df = load_mmpbsa_dat('test/results/6t1h_1151_comp/analysis_20260207_223534/MMPBSA.py_results/FINAL_DECOMP_MMPBSA.dat')
    
    # Apply Mapping: MMPBSA Residue X = OpenGBSA Residue (X + 33)
    # Check if this holds for SER 1 -> SER 34
    
    ser34 = og_res_df[og_res_df['residue_number'] == 34]
    if not ser34.empty:
        print("Mapping Check: Found OpenGBSA Residue 34:", ser34.iloc[0]['residue_name'])
        if ser34.iloc[0]['residue_name'] == 'SER':
            print("Mapping Verified: MMPBSA Res 1 (SER) aligns with OpenGBSA Res 34 (SER).")
            offset = 33
        else:
            print("Mapping Mismatch: Res 34 is", ser34.iloc[0]['residue_name'])
            offset = 0 # Fallback or need investigation
    else:
        print("Mapping Warning: Residue 34 not found in OpenGBSA results.")
        offset = 33 # Assume alignment based on PDB view earlier
        
    mmpbsa_res_df['opengbsa_res_num'] = mmpbsa_res_df['residue_number'] + offset
    
    # Merge
    merged = pd.merge(mmpbsa_res_df, og_res_df, left_on='opengbsa_res_num', right_on='residue_number', how='inner')
    
    # Note: OpenGBSA `total` in per_residue_detailed.csv might be `complex_total` - `receptor_total` - `ligand_total`? 
    # Or is `total` column the Delta?
    # Usually `total` in the CSV is the Delta Total Contribution.
    # Let's verify with values.
    
    # Calculate Differences
    merged['diff_total'] = merged['total'] - merged['mmpbsa_total']
    merged['diff_vdw'] = merged['vdw'] - merged['mmpbsa_vdw']
    merged['diff_ele'] = merged['electrostatic'] - merged['mmpbsa_ele']
    merged['diff_solv'] = merged['solvation'] - (merged['mmpbsa_gb'] + merged['mmpbsa_sa'])
    
    # Sort by absolute diff total
    merged['abs_diff_total'] = merged['diff_total'].abs()
    top_diffs = merged.sort_values('abs_diff_total', ascending=False).head(10)
    
    print(f"\n{'Residue':<15} {'OG_Total':<10} {'MM_Total':<10} {'Diff':<10} {'OG_Solv':<10} {'MM_Solv(GB+SA)':<15}")
    for _, row in top_diffs.iterrows():
        res_label = f"{row['residue_name_x']}_{row['residue_number_x']}({row['residue_number_y']})"
        mm_solv = row['mmpbsa_gb'] + row['mmpbsa_sa']
        print(f"{res_label:<15} {row['total']:<10.2f} {row['mmpbsa_total']:<10.2f} {row['diff_total']:<10.2f} {row['solvation']:<10.2f} {mm_solv:<15.2f}")

    # Check Solvation Zero
    zero_solv_count = (og_res_df['solvation'] == 0.0).sum()
    print(f"\nOpenGBSA Residues with 0.0 Solvation: {zero_solv_count}/{len(og_res_df)}")

except Exception as e:
    print(f"Error analyzing per-residue results: {e}")


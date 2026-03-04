
import pandas as pd
import numpy as np
import os

def check_entropy():
    # Load the results CSV
    csv_file = "test/results/full_feature_test/analysis_20260207_203737/fixed_enhanced_mmgbsa_results_obc2.csv"
    
    if not os.path.exists(csv_file):
        print(f"File not found: {csv_file}")
        return

    df = pd.read_csv(csv_file)
    print(f"Loaded {len(df)} frames.")
    
    if 'total' in df.columns:
        print("Using 'total' column for Binding Energy.")
        delta_H = df['total'].values
    elif 'binding_energy' in df.columns:
        print("Using 'binding_energy' column.")
        delta_H = df['binding_energy'].values
    else:
        print("Could not find binding energy column. Columns:", df.columns)
        return

    # Constants
    temp = 300.0
    R = 0.001987
    RT = R * temp
    beta = 1.0 / RT
    
    # 1. Statistics
    mean_H = np.mean(delta_H)
    std_H = np.std(delta_H, ddof=1)
    print(f"Mean Binding Energy (dH): {mean_H:.2f} kcal/mol")
    print(f"Std Dev (sigma): {std_H:.2f} kcal/mol")
    
    # 2. Gaussian Entropy
    gaussian_penalty = (std_H ** 2) / (2 * RT)
    print(f"Gaussian Entropy Penalty (-TdS): {gaussian_penalty:.2f} kcal/mol")
    print(f"Est. Free Energy (Gaussian): {mean_H + gaussian_penalty:.2f} kcal/mol")
    
    # 3. Interaction Entropy
    try:
        dev = delta_H - mean_H
        arg = beta * dev
        max_arg = np.max(arg)
        sum_exp = np.sum(np.exp(arg - max_arg))
        log_sum = max_arg + np.log(sum_exp)
        log_mean = log_sum - np.log(len(arg))
        ie_val = RT * log_mean
        print(f"Interaction Entropy Penalty (-TdS): {ie_val:.2f} kcal/mol")
        print(f"Est. Free Energy (IE): {mean_H + ie_val:.2f} kcal/mol")
    except Exception as e:
        print(f"IE Calculation Failed: {e}")

if __name__ == "__main__":
    check_entropy()

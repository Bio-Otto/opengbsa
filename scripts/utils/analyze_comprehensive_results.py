
import os
import pandas as pd
import glob
import sys

def analyze_results():
    base_dir = "test/results/comprehensive"
    results = []

    # specific tests to compare for GB model differences
    gb_tests = ['test_gb_obc1', 'test_gb_obc2', 'test_gb_gbn', 'test_gb_gbn2', 'test_gb_hct']
    
    print(f"{'Test Name':<25} | {'Binding Energy':<15} | {'Delta GB':<15} | {'Delta SA':<15}")
    print("-" * 80)

    for test_name in gb_tests:
        # Find the result CSV
        search_path = os.path.join(base_dir, test_name, "analysis_*", "fixed_enhanced_mmgbsa_results_*.csv")
        files = glob.glob(search_path)
        
        if not files:
            print(f"{test_name:<25} | {'MISSING':<15} | {'N/A':<15} | {'N/A':<15}")
            continue
            
        # Sort by mtime to get latest (though there should only be one relevant run)
        latest_file = max(files, key=os.path.getmtime)
        
        try:
            df = pd.read_csv(latest_file)
            # Get mean binding energy
            if 'binding_energy' in df.columns:
                binding = df['binding_energy'].mean()
                gb = df['delta_gb'].mean()
                sa = df['delta_sa'].mean()
                
                print(f"{test_name:<25} | {binding:<15.4f} | {gb:<15.4f} | {sa:<15.4f}")
                results.append(binding)
            else:
                print(f"{test_name:<25} | CSV format error")
                
        except Exception as e:
            print(f"{test_name:<25} | Error reading CSV: {e}")

    # Check for identical values (which would indicate failure)
    print("\n--- Uniqueness Check ---")
    if len(results) < len(gb_tests):
        print(f"WARNING: Only analyzed {len(results)}/{len(gb_tests)} tests. Check missing outputs.")
        sys.exit(1)

    unique_energies = set([round(x, 4) for x in results])
    if len(unique_energies) == len(results):
        print("SUCCESS: All GB models produced distinct binding energies.")
    else:
        print(f"FAILURE: Found duplicate energies! ({len(results)} tests, {len(unique_energies)} unique values)")
        print(results)

if __name__ == "__main__":
    analyze_results()

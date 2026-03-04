
import os
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import shutil
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import sys
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import yaml
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pandas
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
 as pd
import numpy
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
 as np

from mmgbsa.runner import MMGBSARunner
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
 as Runner

# Define paths
# Define paths
DATA_DIR = 'test/data/7khz_gro_test'
TOPOLOGY = os.path.join(DATA_DIR, 'topol_dry.prmtop')
TRAJECTORY = os.path.join(DATA_DIR, 'md_complex_full.xtc')
SOLVATED_TOP = os.path.join(DATA_DIR, 'topol_solvated.prmtop')
LIGAND_RES = 'LIG'

def run_analysis(cutoff, output_dir):
    print(f"Running analysis with nonbonded_cutoff = {cutoff} A...")
    
    config = {
        'input_files': {
            'complex_pdb': TOPOLOGY,
            'trajectory': TRAJECTORY,
            'solvated_topology': SOLVATED_TOP,
            'ligand_resname': LIGAND_RES
        },
        'output_settings': {
            'output_directory': output_dir,
            'overwrite': True,
            'verbose': False
        },
        'analysis_settings': {
            'gb_model': 'OBC2',
            'salt_concentration': 0.15,
            'nonbonded_cutoff': cutoff,
            'run_per_residue_decomposition': True,
            'decomposition_method': 'per_residue',
            'start_frame': 0,
            'end_frame': 2,
            'frame_stride': 1
        }
    }
    
    runner = Runner(config, output_dir=output_dir)
    runner.run_analysis()
    
    # Load decomposition results
    # file matches verify_cutoff_8/frame_by_frame_decomposition.csv
    decomp_file = os.path.join(output_dir, 'frame_by_frame_decomposition.csv')
    
    if not os.path.exists(decomp_file):
        # Fallback to main results if decomp not found
        gbsa_file = os.path.join(output_dir, 'fixed_enhanced_mmgbsa_results_obc2.csv')
        if os.path.exists(gbsa_file):
            print(f"WARNING: Decomposition file not found, using GBSA results: {gbsa_file}")
            return pd.read_csv(gbsa_file)
        
        # List dir for debugging
        print(f"ERROR: Files in {output_dir}: {os.listdir(output_dir)}")
        raise FileNotFoundError(f"Decomposition file not found: {decomp_file}")
        
    df = pd.read_csv(decomp_file)
    return df

def main():
    print("Verifying Decomposition Cutoff Sensitivity with OpenMM 8.5...")
    
    out_8 = 'verify_cutoff_8'
    out_12 = 'verify_cutoff_12'
    
    try:
        df_8 = run_analysis(8.0, out_8)
        df_12 = run_analysis(12.0, out_12)
        
        # Compare Total Energy
        # frame_by_frame_decomposition.csv has columns like "RES_Total", etc?
        # Or just Residues. 
        # Let's sum all numeric columns to get "Total" if "Total" column missing
        
        if 'Total' in df_8.columns:
            total_8 = df_8['Total'].values
            total_12 = df_12['Total'].values
        else:
            # Sum all residue columns (assuming frame index is not a column or handled)
            # actually usually valid columns are numeric.
            print("Calculating Total from residue sum...")
            total_8 = df_8.select_dtypes(include=[np.number]).sum(axis=1).values
            total_12 = df_12.select_dtypes(include=[np.number]).sum(axis=1).values

        diff = np.abs(total_8 - total_12)
        max_diff = np.max(diff)
        mean_diff = np.mean(diff)
        
        print(f"\nComparing Results:")
        print(f"Max Difference in Total Energy: {max_diff:.4f} kcal/mol")
        print(f"Mean Difference in Total Energy: {mean_diff:.4f} kcal/mol")
        
        if max_diff > 0.01:
            print("\nSUCCESS: Decomposition results match expectations (Sensitivity detected).")
            print("The fix for nonbonded_cutoff is active and working.")
        else:
            print("\nFAILURE: Decomposition results are identical!")
            print("The nonbonded_cutoff parameter is NOT affecting the results.")
            sys.exit(1)
            
    except Exception as e:
        print(f"\nERROR: {e}")
        import traceback
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

        traceback.print_exc()
        sys.exit(1)
    finally:
        # Cleanup - DISABLED FOR DEBUGGING
        pass
        # if os.path.exists(out_8):
        #     shutil.rmtree(out_8)
        # if os.path.exists(out_12):
        #     shutil.rmtree(out_12)

if __name__ == "__main__":
    main()

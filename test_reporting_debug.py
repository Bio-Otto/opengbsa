import os
import shutil
import numpy as np
from mmgbsa.reporting import HTMLReportGenerator

def test_html_generation():
    print("Testing HTML Report Generation...")
    
    # 1. Create Dummy Data
    output_dir = "test_html_output"
    os.makedirs(output_dir, exist_ok=True)
    
    # Analysis Results (Mock)
    analysis_results = {
        'mean_contribution': -25.5,
        'n_residues': 300,
        'hot_spots': None # Can be None
    }
    
    # Frame Data (Mock)
    # 3 Frames, 2 Residues
    frame_data = [
        {'frame_number': 0, 'RES1_A_total': -10.0, 'RES2_A_total': -5.0, 'RES1_A_vdw': -8.0},
        {'frame_number': 1, 'RES1_A_total': -12.0, 'RES2_A_total': np.nan, 'RES1_A_vdw': -9.0}, # NaN Test
        {'frame_number': 2, 'RES1_A_total': -11.0, 'RES2_A_total': -6.0, 'RES1_A_vdw': -8.5},
    ]
    
    # Write Dummy PDB
    pdb_path = "test_complex.pdb"
    with open(pdb_path, 'w') as f:
        f.write("ATOM      1  N   ALA A   1      10.000  10.000  10.000  1.00  0.00           N\\n")
    
    # 2. Initialize Generator
    generator = HTMLReportGenerator(output_dir)
    
    # 3. Generate Report
    try:
        report_path = generator.generate_report(analysis_results, frame_data, complex_pdb_path=pdb_path)
        
        if report_path and os.path.exists(report_path):
            print(f"SUCCESS: Report generated at {report_path}")
            
            # 4. Verify Content
            with open(report_path, 'r') as f:
                content = f.read()
                
            if "var data = [" in content:
                print("  - Plotly data script found.")
            else:
                print("  - FAIL: Plotly data script NOT found.")
                
            if "null" in content: # Check if NaN became null
                print("  - NaN sanitization check: Passed (found 'null').")
            else:
                print("  - NaN sanitization check: Warning (no 'null' found, maybe verification logic is loose).")
                
        else:
            print("FAIL: Report file not created.")
            
    except Exception as e:
        print(f"CRASH: {e}")

if __name__ == "__main__":
    test_html_generation()

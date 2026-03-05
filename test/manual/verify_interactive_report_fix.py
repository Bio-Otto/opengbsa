
import os
import pandas as pd
import numpy as np
from mmgbsa.reporting import HTMLReportGenerator

def test_report_generation_with_config():
    print("Testing HTML Report Generation with Entropy Config...")
    
    # Mock Data
    frame_data = [{'frame_index': 0, 'frame_number': 1, 'ALA1_total': -1.5}]
    
    # Mock Global Results (Binding Energy series for entropy calc)
    # Create variance to ensure entropy term is non-zero
    global_df = pd.DataFrame({
        'frame': range(1, 11),
        'binding_energy': np.random.normal(-30, 2, 10), # Mean -30, Std 2
        'delta_nb': np.random.normal(-40, 2, 10),
        'delta_gb': np.random.normal(15, 1, 10),
        'delta_sa': np.random.normal(-5, 0.1, 10),
        'delta_vdw': np.random.normal(-35, 1, 10),
        'delta_elec': np.random.normal(-5, 1, 10)
    })
    
    global_results = {'dataframe': global_df}
    analysis_results = {'total_contribution': -30.0}
    
    # Mock Config with Analysis Settings (testing fallback)
    config = {
        'analysis_settings': {
            'entropy_method': 'interaction'
        },
        # reporting_settings intentionally left empty or without entropy_approximation
    }
    
    output_dir = "test/results/report_verification_entropy"
    os.makedirs(output_dir, exist_ok=True)
    
    generator = HTMLReportGenerator(output_dir, config=config)
    
    try:
        # Generate
        report_path = generator.generate_report(
            analysis_results, 
            frame_data, 
            global_results=global_results, 
            complex_pdb_path="dummy.pdb", 
            ligand_resname="LIG"
        )
        
        if report_path and os.path.exists(report_path):
            print(f"✅ Report generated: {report_path}")
            
            with open(report_path, 'r') as f:
                content = f.read()
                
                # Check 1: Does it mention "Interaction Entropy"?
                if "Interaction Entropy" in content:
                    print("✅ Found 'Interaction Entropy' in report (Config propagation successful).")
                else:
                    print("❌ 'Interaction Entropy' NOT found. Fallback logic failed.")
                    
                # Check 2: Does "Total Free Energy (ΔG)" appear in the Bar Chart section?
                # The label "Total Free Energy (ΔG)" should produce a trace name or axis label.
                if "Total Free Energy (ΔG)" in content:
                    print("✅ Found 'Total Free Energy (ΔG)' label (Bar chart update successful).")
                else:
                    print("❌ 'Total Free Energy (ΔG)' label NOT found.")
                    
        else:
             print("❌ Report generation failed.")
             
    except Exception as e:
        print(f"❌ Exception: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    test_report_generation_with_config()

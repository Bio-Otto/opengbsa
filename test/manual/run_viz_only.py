
import pandas as pd
import numpy as np
from pathlib import Path
from mmgbsa.visualization import AdvancedVisualization
import warnings
warnings.filterwarnings('ignore')

# Path to existing results
output_dir = Path("test/results/full_feature_test/analysis_20260207_174017")
print(f"Loading results from {output_dir}")

# Load Data
binding_df = pd.read_csv(output_dir / "fixed_enhanced_mmgbsa_results_obc2.csv")
decomp_df = pd.read_csv(output_dir / "per_residue_detailed.csv")

# Ensure proper types
if 'frame' not in binding_df.columns:
    binding_df['frame'] = range(len(binding_df)) 

# Construct input matching core.py logic
# Note: we use delta_g (mean binding energy) for IC50
mean_dg = binding_df['binding_energy'].mean()
std_dg = binding_df['binding_energy'].std()


print(f"Mean Binding Energy: {mean_dg:.2f} +/- {std_dg:.2f}")

input_data = {
    'mean_binding_energy': mean_dg,
    'std_binding_energy': std_dg,
    'temperature': 300.0,
    'n_frames': len(binding_df),
    'binding_data': binding_df,
    'decomposition_results': {'dataframe': decomp_df}
}

# Run Visualization
viz = AdvancedVisualization(output_dir)
viz.load_mmgbsa_results(input_data)
viz.generate_comprehensive_plots(compound_name="Ligand_7khz")
viz.save_results()

# Also generate the interactive report from reporting.py
from mmgbsa.reporting import HTMLReportGenerator
print("\nAlso generating standard interactive report...")
output_dir = Path("test/results/full_feature_test/analysis_20260207_174017")
input_data['total_contribution'] = -32.61 # Manually set for the test
input_data['dataframe'] = binding_df # Ensure binding data is available
report_gen = HTMLReportGenerator(str(output_dir), config={ 'reporting_settings': {'entropy_approximation': 'gaussian'}})
report_gen.generate_report(input_data, binding_df.to_dict('records'), global_results=input_data)

print("\nVisualization complete.")

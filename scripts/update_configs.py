import os
import glob
import yaml
import random

configs = glob.glob('test/configs/comprehensive/*.yaml')

updated_count = 0
for config_file in configs:
    with open(config_file, 'r') as f:
        try:
            data = yaml.safe_load(f)
        except Exception as e:
            print(f"Error loading {config_file}: {e}")
            continue
    
    if data is None: 
        continue

    # Ensure sections exist
    if 'analysis_settings' not in data:
        data['analysis_settings'] = {}
    if 'platform_settings' not in data:
        data['platform_settings'] = {}

    # 1. Decomposition settings (at least 10 frames)
    data['analysis_settings']['run_per_residue_decomposition'] = True
    
    current_decomp = data['analysis_settings'].get('decomp_frames', 0)
    if current_decomp is None:
        current_decomp = 0
    new_decomp = max(10, current_decomp)
    data['analysis_settings']['decomp_frames'] = new_decomp

    # Ensure max_frames >= decomp_frames
    current_max = data['analysis_settings'].get('max_frames', 0)
    if current_max is None:
        current_max = 0
    data['analysis_settings']['max_frames'] = max(current_max, new_decomp)

    # 2. CPU multiprocess decompose
    data['analysis_settings']['parallel_processing'] = True
    data['platform_settings']['decomposition_platform'] = 'CPU'

    # 3. Randomize CUDA/OpenCL for main
    data['platform_settings']['preferred_platform'] = random.choice(['CUDA', 'OpenCL'])

    # Write the updated config back
    with open(config_file, 'w') as f:
        yaml.dump(data, f, default_flow_style=False, sort_keys=False)
    
    updated_count += 1

print(f"Successfully updated {updated_count} configuration files.")

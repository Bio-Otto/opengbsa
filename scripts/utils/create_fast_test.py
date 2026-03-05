import yaml

with open('6xj3_analysis_config.yaml', 'r') as f:
    config = yaml.safe_load(f)

# Change max frames to 5 for speed
config['analysis_settings']['max_frames'] = 5
config['analysis_settings']['decomp_frames'] = 5
config['analysis_settings']['frame_start'] = 1000
config['analysis_settings']['frame_stride'] = 1

# Output to a specific test dir
config['output_settings']['output_directory'] = 'test_fix_output'

with open('fast_test_config.yaml', 'w') as f:
    yaml.dump(config, f)

print("Created fast_test_config.yaml. Now you can run it.")

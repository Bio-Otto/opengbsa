
import os
import yaml
from copy import deepcopy

# Base Configuration (Template)
BASE_CONFIG = {
    'input_files': {
        'complex_pdb': 'test/data/6t1h_prmtop_test/complex.prmtop',
        'trajectory': 'test/data/6t1h_prmtop_test/md_complex_prod.xtc',
        'solvated_topology': 'test/data/6t1h_prmtop_test/replica_1/gro/complex_solv_ions.gro',
        'receptor_topology': 'test/data/6t1h_prmtop_test/receptor.prmtop',
        'ligand_topology': 'test/data/6t1h_prmtop_test/ligand.prmtop'
    },
    'output_settings': {
        'save_plots': True,
        'save_logs': True
    },
    'analysis_settings': {
        'temperature': 300.0,
        'frame_stride': 1,
        'max_frames': 1,  # Minimal frames for testing (SPEED UP)
        'frame_selection': 'sequential',
        'use_cache': False, # Disable cache to force calculation
        'verbose': 1,
        'run_entropy_analysis': False,
        'entropy_method': 'none',
        'run_per_residue_decomposition': False,
        'report_raw_energies': False,
        'decomp_frames': 1  # Minimal frames if decomposition enabled
    },
    'forcefield_settings': {
        'protein_forcefield': 'amber'
    }
}

OUTPUT_DIR = 'test/configs/comprehensive'
os.makedirs(OUTPUT_DIR, exist_ok=True)

test_cases = []

# 1. Test All GB Models (Base settings: ACE, Salt 0.15)
gb_models = ['OBC1', 'OBC2', 'GBn', 'GBn2', 'HCT']
for gb in gb_models:
    test_cases.append({
        'name': f'test_gb_{gb.lower()}',
        'settings': {'gb_model': gb, 'salt_concentration': 0.15, 'sa_model': 'ACE'}
    })

# 2. Test Salt Concentrations with GBn
salts = [0.0, 0.5]
for salt in salts:
    test_cases.append({
        'name': f'test_salt_{str(salt).replace(".", "")}',
        'settings': {'gb_model': 'GBn', 'salt_concentration': salt, 'sa_model': 'ACE'}
    })

# 3. Test SA Model LCPO with OBC2
test_cases.append({
    'name': 'test_sa_lcpo',
    'settings': {'gb_model': 'OBC2', 'salt_concentration': 0.15, 'sa_model': 'LCPO'}
})

# 4. Test Entropy Methods with GBn
entropies = ['interaction', 'quasiharmonic'] # Skipping normal_mode for speed unless requested
for ent in entropies:
    test_cases.append({
        'name': f'test_entropy_{ent}',
        'settings': {
            'gb_model': 'GBn', 
            'run_entropy_analysis': True, 
            'entropy_method': ent
        }
    })

# 5. Test Decomposition with Report Raw Energies
test_cases.append({
    'name': 'test_decomp_raw',
    'settings': {
        'gb_model': 'GBn',
        'run_per_residue_decomposition': True,
        'report_raw_energies': True
    }
})

# 6. Test Decomposition without Raw Energies
test_cases.append({
    'name': 'test_decomp_noraw',
    'settings': {
        'gb_model': 'OBC2',
        'run_per_residue_decomposition': True,
        'report_raw_energies': False
    }
})

# 7. Additional Combinations (Mix and Match)
combinations = [
    # HCT + Interaction Entropy
    {'name': 'test_hct_ie', 'settings': {'gb_model': 'HCT', 'run_entropy_analysis': True, 'entropy_method': 'interaction'}},
    # OBC1 + Quasiharmonic
    {'name': 'test_obc1_qha', 'settings': {'gb_model': 'OBC1', 'run_entropy_analysis': True, 'entropy_method': 'quasiharmonic'}},
    # GBn2 + Salt 0.0 + Decomp
    {'name': 'test_gbn2_nosalt_decomp', 'settings': {'gb_model': 'GBn2', 'salt_concentration': 0.0, 'run_per_residue_decomposition': True, 'report_raw_energies': True}},
    # OBC2 + LCPO + Decomp
    {'name': 'test_obc2_lcpo_decomp', 'settings': {'gb_model': 'OBC2', 'sa_model': 'LCPO', 'run_per_residue_decomposition': True, 'report_raw_energies': True}},
]
test_cases.extend(combinations)


# Generate Files
print(f"Generating {len(test_cases)} test configurations...")

for case in test_cases:
    config = deepcopy(BASE_CONFIG)
    
    # Update Analysis Name and Output Directory
    config['output_settings']['analysis_name'] = case['name']
    config['output_settings']['output_directory'] = f'test/results/comprehensive/{case["name"]}'
    
    # Update specific settings
    # Default settings that might be missing in case['settings'] but needed default if not in BASE
    config['analysis_settings']['gb_model'] = 'GBn' # Default
    config['analysis_settings']['salt_concentration'] = 0.15 # Default
    config['analysis_settings']['sa_model'] = 'ACE' # Default

    config['analysis_settings'].update(case['settings'])
    
    # Write to file
    filename = f"{OUTPUT_DIR}/{case['name']}.yaml"
    with open(filename, 'w') as f:
        yaml.dump(config, f, default_flow_style=False)
    
    print(f"Created: {filename}")

print("Done.")

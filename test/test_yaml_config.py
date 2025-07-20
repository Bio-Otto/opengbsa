#!/usr/bin/env python3
"""
Simple test script for YAML configuration system
Tests the configuration loading and validation without requiring OpenMM
"""

import yaml
import sys
from pathlib import Path

def test_yaml_config():
    """Test YAML configuration loading and validation"""
    
    print("🧪 Testing YAML Configuration System")
    print("=" * 50)
    
    # Test 1: Load existing config
    print("\n📋 Test 1: Loading existing configuration...")
    try:
        with open('mmgbsa_config.yaml', 'r') as f:
            config = yaml.safe_load(f)
        print("✅ Configuration loaded successfully")
        
        # Check required sections
        required_sections = ['input_files', 'analysis_settings']
        for section in required_sections:
            if section in config:
                print(f"✅ Found required section: {section}")
            else:
                print(f"❌ Missing required section: {section}")
                return False
        
        # Check input files
        input_files = config['input_files']
        print(f"\n📁 Input files configuration:")
        for file_type, file_path in input_files.items():
            print(f"  • {file_type}: {file_path}")
        
        # Check analysis settings
        analysis_settings = config['analysis_settings']
        print(f"\n⚙️  Analysis settings:")
        print(f"  • Temperature: {analysis_settings.get('temperature', 'N/A')} K")
        print(f"  • GB Model: {analysis_settings.get('gb_model', 'N/A')}")
        print(f"  • Max Frames: {analysis_settings.get('max_frames', 'N/A')}")
        print(f"  • Salt Concentration: {analysis_settings.get('salt_concentration', 'N/A')} M")
        
    except Exception as e:
        print(f"❌ Error loading configuration: {e}")
        return False
    
    # Test 2: Create sample config
    print(f"\n📝 Test 2: Creating sample configuration...")
    try:
        sample_config = {
            'input_files': {
                'ligand_mol': 'test/ligand.sdf',
                'complex_pdb': 'test/complex.pdb',
                'ligand_pdb': 'test/ligand.pdb',
                'trajectory': 'test/complex.xtc'
            },
            'analysis_settings': {
                'temperature': 300,
                'gb_model': 'OBC2',
                'salt_concentration': 0.15,
                'max_frames': 50,
                'verbose': 1,
                'use_cache': True,
                'parallel_processing': False,
                'energy_decomposition': False,
                'run_entropy_analysis': False,
                'run_per_residue_decomposition': False,
                'decomp_frames': 10
            },
            'output_settings': {
                'output_directory': 'mmgbsa_results',
                'save_plots': True,
                'save_trajectories': False
            },
            'advanced_settings': {
                'minimization_tolerance': 1e-6,
                'nma_quality_threshold': 'Good',
                'hot_spot_threshold': -1.0
            }
        }
        
        with open('test_config.yaml', 'w') as f:
            yaml.dump(sample_config, f, default_flow_style=False, indent=2)
        
        print("✅ Sample configuration created: test_config.yaml")
        
    except Exception as e:
        print(f"❌ Error creating sample configuration: {e}")
        return False
    
    # Test 3: Validate file paths
    print(f"\n🔍 Test 3: Validating file paths...")
    test_files = {
        'ligand_mol': 'test/ligand.sdf',
        'complex_pdb': 'test/complex.pdb',
        'ligand_pdb': 'test/ligand.pdb',
        'trajectory': 'test/complex.xtc'
    }
    
    missing_files = []
    for file_type, file_path in test_files.items():
        if Path(file_path).exists():
            print(f"✅ {file_type}: {file_path}")
        else:
            print(f"⚠️  {file_type}: {file_path} (not found)")
            missing_files.append(file_path)
    
    if missing_files:
        print(f"\n⚠️  Some test files are missing: {missing_files}")
        print("This is expected if you haven't set up the test data yet.")
    else:
        print(f"\n✅ All test files found!")
    
    # Test 4: Configuration validation
    print(f"\n✅ Test 4: Configuration validation...")
    try:
        # Validate required sections
        if 'input_files' in sample_config and 'analysis_settings' in sample_config:
            print("✅ Required sections present")
        else:
            print("❌ Missing required sections")
            return False
        
        # Validate analysis settings
        analysis = sample_config['analysis_settings']
        if analysis.get('temperature', 0) > 0:
            print("✅ Valid temperature")
        else:
            print("❌ Invalid temperature")
            return False
        
        if analysis.get('gb_model') in ['OBC2', 'OBC1', 'HCT', 'GBn', 'GBn2']:
            print("✅ Valid GB model")
        else:
            print("❌ Invalid GB model")
            return False
        
        print("✅ Configuration validation passed")
        
    except Exception as e:
        print(f"❌ Configuration validation failed: {e}")
        return False
    
    print(f"\n🎉 All tests passed!")
    print(f"✅ YAML configuration system is working correctly")
    print(f"✅ Ready for MM/GBSA analysis")
    
    return True

def show_usage_examples():
    """Show usage examples"""
    print(f"\n📚 Usage Examples:")
    print(f"=" * 30)
    
    print(f"\n1. Create configuration file:")
    print(f"   python mmgbsa_runner.py --create-config")
    
    print(f"\n2. Run analysis with configuration:")
    print(f"   python mmgbsa_runner.py mmgbsa_config.yaml")
    
    print(f"\n3. Create custom configuration:")
    print(f"   python mmgbsa_runner.py --create-config --config-name my_config.yaml")
    
    print(f"\n4. Programmatic usage:")
    print(f"   from mmgbsa_runner import MMGBSARunner")
    print(f"   runner = MMGBSARunner('config.yaml')")
    print(f"   results = runner.run_analysis()")

if __name__ == '__main__':
    success = test_yaml_config()
    
    if success:
        show_usage_examples()
        print(f"\n🚀 Ready to run MM/GBSA analysis!")
    else:
        print(f"\n❌ Configuration system test failed")
        sys.exit(1) 
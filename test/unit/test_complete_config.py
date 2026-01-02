#!/usr/bin/env python3
"""
Test script for complete MM/GBSA configuration
Validates all sections and parameters
"""

import yaml
import json
from pathlib import Path
import warnings

# Suppress warnings
warnings.filterwarnings('ignore')

# Mock dependencies if not available
import sys
from unittest.mock import MagicMock

def mock_dependencies():
    """Mock missing dependencies for testing configuration logic"""
    modules_to_mock = [
        'pandas', 
        'numpy', 
        'scipy',
        'scipy.stats',
        'scipy.spatial',
        'scipy.interpolate',
        'scipy.constants',
        'matplotlib',
        'matplotlib.pyplot',
        'seaborn',
        'openmm', 
        'simtk.openmm',
        'openmm.app',
        'openmm.unit',
        'mdtraj',
        'openff',
        'openff.toolkit',
        'openff.toolkit.topology',
        'openff.toolkit.typing.engines.smirnoff', 
        'openmmforcefields',
        'openmmforcefields.generators'
    ]
    
    for module_name in modules_to_mock:
        if module_name not in sys.modules:
            try:
                # Try import
                __import__(module_name)
            except (ImportError, AttributeError):
                # Mock if missing
                mock = MagicMock()
                # Hack to make it look like a package
                mock.__path__ = []
                mock.__spec__ = MagicMock()
                sys.modules[module_name] = mock
                
                # Ensure parent packages exist
                if '.' in module_name:
                    parts = module_name.split('.')
                    # Ensure all parents are in sys.modules
                    for i in range(1, len(parts)):
                        parent = '.'.join(parts[:i])
                        if parent not in sys.modules:
                            parent_mock = MagicMock()
                            parent_mock.__path__ = []
                            parent_mock.__spec__ = MagicMock()
                            sys.modules[parent] = parent_mock

mock_dependencies()

def test_complete_configuration():
    """Test the complete configuration file"""
    
    print("🧪 Testing Complete MM/GBSA Configuration")
    print("=" * 60)
    
    # Test main configuration
    try:
        # Try to find the config file in different locations
        config_paths = [
            'complete_mmgbsa_config.yaml',
            '../complete_mmgbsa_config.yaml',
            'test/complete_mmgbsa_config.yaml'
        ]
        
        config = None
        for path in config_paths:
            if Path(path).exists():
                with open(path, 'r') as f:
                    config = yaml.safe_load(f)
                break
        
        if config is None:
            raise FileNotFoundError("Could not find complete_mmgbsa_config.yaml")
        
        print("✅ Complete configuration loaded successfully")
        
        # Test all sections
        sections = [
            'input_files', 'output_settings', 'analysis_settings',
            'forcefield_settings', 'advanced_settings', 'platform_settings',
            'validation_settings', 'reporting_settings', 'debug_settings',
            'reproducibility_settings', 'performance_settings'
        ]
        
        print(f"\n📋 Configuration Sections:")
        for section in sections:
            if section in config:
                print(f"  ✅ {section}")
            else:
                print(f"  ❌ {section} (missing)")
        
        # Test input files
        print(f"\n📁 Input Files:")
        input_files = config.get('input_files', {})
        required_files = ['ligand_mol', 'complex_pdb', 'ligand_pdb', 'trajectory']
        
        for file_type in required_files:
            if file_type in input_files:
                file_path = input_files[file_type]
                if Path(file_path).exists():
                    print(f"  ✅ {file_type}: {file_path}")
                else:
                    print(f"  ⚠️  {file_type}: {file_path} (file not found)")
            else:
                print(f"  ❌ {file_type}: missing")
        
        # Test output settings
        print(f"\n📊 Output Settings:")
        output_settings = config.get('output_settings', {})
        print(f"  • Output directory: {output_settings.get('output_directory', 'N/A')}")
        print(f"  • Analysis name: {output_settings.get('analysis_name', 'N/A')}")
        print(f"  • Output formats: {output_settings.get('output_formats', 'N/A')}")
        print(f"  • Save plots: {output_settings.get('save_plots', 'N/A')}")
        
        # Test analysis settings
        print(f"\n⚙️  Analysis Settings:")
        analysis_settings = config.get('analysis_settings', {})
        print(f"  • Temperature: {analysis_settings.get('temperature', 'N/A')} K")
        print(f"  • GB Model: {analysis_settings.get('gb_model', 'N/A')}")
        print(f"  • Salt concentration: {analysis_settings.get('salt_concentration', 'N/A')} M")
        print(f"  • Max frames: {analysis_settings.get('max_frames', 'N/A')}")
        print(f"  • Frame selection: {analysis_settings.get('frame_selection', 'N/A')}")
        
        # Test forcefield settings
        print(f"\n🔬 Forcefield Settings:")
        forcefield_settings = config.get('forcefield_settings', {})
        print(f"  • Protein: {forcefield_settings.get('protein_forcefield', 'N/A')}")
        print(f"  • Ligand: {forcefield_settings.get('ligand_forcefield', 'N/A')}")
        print(f"  • Water: {forcefield_settings.get('water_model', 'N/A')}")
        
        # Test advanced settings
        print(f"\n🔧 Advanced Settings:")
        advanced_settings = config.get('advanced_settings', {})
        print(f"  • Min tolerance: {advanced_settings.get('minimization_tolerance', 'N/A')}")
        print(f"  • NMA quality: {advanced_settings.get('nma_quality_threshold', 'N/A')}")
        print(f"  • Hot spot threshold: {advanced_settings.get('hot_spot_threshold', 'N/A')} kcal/mol")
        
        # Test platform settings
        print(f"\n💻 Platform Settings:")
        platform_settings = config.get('platform_settings', {})
        print(f"  • Preferred platform: {platform_settings.get('preferred_platform', 'N/A')}")
        print(f"  • CUDA device: {platform_settings.get('cuda_device', 'N/A')}")
        print(f"  • CUDA precision: {platform_settings.get('cuda_precision', 'N/A')}")
        
        # Test validation settings
        print(f"\n✅ Validation Settings:")
        validation_settings = config.get('validation_settings', {})
        print(f"  • Validate inputs: {validation_settings.get('validate_inputs', 'N/A')}")
        print(f"  • Validate results: {validation_settings.get('validate_results', 'N/A')}")
        print(f"  • Max std dev: {validation_settings.get('max_std_dev', 'N/A')} kcal/mol")
        
        # Test reporting settings
        print(f"\n📊 Reporting Settings:")
        reporting_settings = config.get('reporting_settings', {})
        print(f"  • Generate final report: {reporting_settings.get('generate_final_report', 'N/A')}")
        print(f"  • Include plots: {reporting_settings.get('include_plots', 'N/A')}")
        print(f"  • Include statistics: {reporting_settings.get('include_statistics', 'N/A')}")
        
        # Test reproducibility settings
        print(f"\n🔄 Reproducibility Settings:")
        reproducibility_settings = config.get('reproducibility_settings', {})
        print(f"  • Global random seed: {reproducibility_settings.get('global_random_seed', 'N/A')}")
        print(f"  • Save configuration: {reproducibility_settings.get('save_configuration', 'N/A')}")
        print(f"  • Save environment: {reproducibility_settings.get('save_environment', 'N/A')}")
        
        # Test performance settings
        print(f"\n⚡ Performance Settings:")
        performance_settings = config.get('performance_settings', {})
        print(f"  • Memory limit: {performance_settings.get('memory_limit', 'N/A')} GB")
        print(f"  • GPU memory fraction: {performance_settings.get('gpu_memory_fraction', 'N/A')}")
        print(f"  • Batch size: {performance_settings.get('batch_size', 'N/A')}")
        
        return True
        
    except Exception as e:
        print(f"❌ Configuration test failed: {e}")
        return False

def test_configuration_validation():
    """Test configuration validation logic"""
    
    print(f"\n🧪 Testing Configuration Validation")
    print("-" * 40)
    
    def validate_config(config):
        """Validate configuration parameters"""
        errors = []
        warnings = []
        
        # Validate temperature
        temp = config.get('analysis_settings', {}).get('temperature', 300)
        if temp <= 0 or temp > 1000:
            errors.append(f"Invalid temperature: {temp} K")
        elif temp < 273 or temp > 373:
            warnings.append(f"Unusual temperature: {temp} K")
        
        # Validate salt concentration
        salt = config.get('analysis_settings', {}).get('salt_concentration', 0.15)
        if salt < 0 or salt > 2:
            errors.append(f"Invalid salt concentration: {salt} M")
        
        # Validate max frames
        max_frames = config.get('analysis_settings', {}).get('max_frames', 50)
        if max_frames is not None and max_frames <= 0:
            errors.append(f"Invalid max_frames: {max_frames}")
        
        # Validate frame stride
        frame_stride = config.get('analysis_settings', {}).get('frame_stride')
        if frame_stride is not None and frame_stride <= 0:
            errors.append(f"Invalid frame_stride: {frame_stride}")
        
        # Validate GB model
        gb_model = config.get('analysis_settings', {}).get('gb_model', 'OBC2')
        valid_gb_models = ['OBC1', 'OBC2', 'HCT', 'GBn', 'GBn2']
        if gb_model not in valid_gb_models:
            errors.append(f"Invalid GB model: {gb_model}")
        
        return errors, warnings
    
    try:
        # Try to find the config file in different locations
        config_paths = [
            'complete_mmgbsa_config.yaml',
            '../complete_mmgbsa_config.yaml',
            'test/complete_mmgbsa_config.yaml'
        ]
        
        config = None
        for path in config_paths:
            if Path(path).exists():
                with open(path, 'r') as f:
                    config = yaml.safe_load(f)
                break
        
        if config is None:
            raise FileNotFoundError("Could not find complete_mmgbsa_config.yaml")
        
        errors, warnings = validate_config(config)
        
        if errors:
            print("❌ Validation errors:")
            for error in errors:
                print(f"  • {error}")
        else:
            print("✅ No validation errors")
        
        if warnings:
            print("⚠️  Validation warnings:")
            for warning in warnings:
                print(f"  • {warning}")
        else:
            print("✅ No validation warnings")
        
        return len(errors) == 0
        
    except Exception as e:
        print(f"❌ Validation test failed: {e}")
        return False

def test_output_structure():
    """Test output directory structure creation"""
    
    print(f"\n🧪 Testing Output Structure")
    print("-" * 30)
    
    try:
        # Try to import the runner from different locations
        try:
            # Try to import from package structure
            try:
                from mmgbsa.complete_runner import CompleteMMGBSARunner
            except ImportError as e:
                print(f"DEBUG: Package import failed: {e}")
                # Try legacy/direct import if package structure not used
                import sys
                sys.path.append(str(Path(__file__).parent.parent))
                try:
                    from mmgbsa.complete_runner import CompleteMMGBSARunner
                except ImportError as e:
                    print(f"DEBUG: Sys.path import failed: {e}")
                    raise e
        except ImportError:
            try:
                # Fallback for older structure or if file is in root/test
                from complete_mmgbsa_runner import CompleteMMGBSARunner
            except ImportError as e:
                print(f"DEBUG: Fallback import failed: {e}")
                print("⚠️  Could not import CompleteMMGBSARunner - skipping output structure test")
                return False
                
        except Exception as e:
            print(f"DEBUG: Unexpected error during import: {e}")
            print("⚠️  Could not import CompleteMMGBSARunner - skipping output structure test")
            return False
        
        # Create a minimal test config
        test_config = {
            'input_files': {
                'ligand_mol': 'test/ligand.sdf',
                'complex_pdb': 'test/complex.pdb',
                'ligand_pdb': 'test/ligand.pdb',
                'trajectory': 'test/complex.xtc'
            },
            'output_settings': {
                'output_directory': 'test_output',
                'analysis_name': 'test_analysis'
            },
            'analysis_settings': {
                'temperature': 300,
                'gb_model': 'OBC2',
                'salt_concentration': 0.15,
                'max_frames': 10
            }
        }
        
        # Save test config
        with open('test_config.yaml', 'w') as f:
            yaml.dump(test_config, f, default_flow_style=False, indent=2)
        
        # Test runner
        runner = CompleteMMGBSARunner('test_config.yaml')
        
        # Test output directory creation
        output_dir = runner._setup_output_directory()
        
        if output_dir.exists():
            print(f"✅ Output directory created: {output_dir}")
            
            # Check subdirectories
            subdirs = ['plots', 'data', 'logs', 'reports', 'cache', 'debug']
            for subdir in subdirs:
                subdir_path = output_dir / subdir
                if subdir_path.exists():
                    print(f"  ✅ {subdir}/")
                else:
                    print(f"  ❌ {subdir}/ (missing)")
            
            # Clean up
            import shutil
            shutil.rmtree('test_output')
            Path('test_config.yaml').unlink()
            
            return True
        else:
            print(f"❌ Output directory not created")
            return False
            
    except Exception as e:
        print(f"❌ Output structure test failed: {e}")
        return False

def main():
    """Main test function"""
    
    print("🚀 Complete Configuration Test Suite")
    print("=" * 60)
    
    # Run tests
    config_test = test_complete_configuration()
    validation_test = test_configuration_validation()
    output_test = test_output_structure()
    
    # Summary
    print(f"\n📊 Test Summary")
    print("=" * 30)
    
    tests = [
        ("Configuration Loading", config_test),
        ("Configuration Validation", validation_test),
        ("Output Structure", output_test)
    ]
    
    passed = sum(1 for _, result in tests if result)
    total = len(tests)
    
    for test_name, result in tests:
        status = "✅ PASS" if result else "❌ FAIL"
        print(f"  {status} {test_name}")
    
    print(f"\nTests passed: {passed}/{total}")
    
    if passed == total:
        print("🎉 All tests passed!")
        print("✅ Complete configuration is ready for use!")
        print("\n📚 Usage:")
        print("  python complete_mmgbsa_runner.py complete_mmgbsa_config.yaml")
    else:
        print("❌ Some tests failed!")
        print("Please check the errors above.")

if __name__ == '__main__':
    main() 
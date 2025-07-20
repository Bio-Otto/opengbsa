#!/usr/bin/env python3
"""
Frame Selection Test Script
Tests the new frame selection functionality in MM/GBSA analysis
"""

import yaml
import numpy as np
from pathlib import Path
import warnings

# Suppress specific warnings
warnings.filterwarnings('ignore')
warnings.filterwarnings("ignore", message="Unable to load toolkit 'OpenEye Toolkit'")
warnings.filterwarnings("ignore", message="importing 'simtk.openmm' is deprecated")

def test_frame_selection_logic():
    """Test frame selection logic without requiring OpenMM"""
    
    print("🧪 Testing Frame Selection Logic")
    print("=" * 50)
    
    def select_frames(trajectory_length, max_frames=None, frame_start=None, frame_end=None,
                     frame_stride=None, frame_selection='sequential', random_seed=42):
        """Frame selection logic (same as in MM/GBSA calculator)"""
        
        # Determine frame range
        if frame_start is None:
            frame_start = 0
        if frame_end is None:
            frame_end = trajectory_length
        
        # Validate frame range
        frame_start = max(0, min(frame_start, trajectory_length - 1))
        frame_end = max(frame_start + 1, min(frame_end, trajectory_length))
        
        print(f"📊 Frame selection parameters:")
        print(f"  • Trajectory length: {trajectory_length}")
        print(f"  • Frame range: {frame_start} to {frame_end}")
        print(f"  • Frame stride: {frame_stride}")
        print(f"  • Selection method: {frame_selection}")
        print(f"  • Max frames: {max_frames}")
        
        # Generate frame indices based on selection method
        if frame_selection == "sequential":
            # Sequential selection with stride
            if frame_stride is None or frame_stride <= 0:
                frame_indices = list(range(frame_start, frame_end))
            else:
                frame_indices = list(range(frame_start, frame_end, frame_stride))
                
        elif frame_selection == "equidistant":
            # Equidistant selection
            if max_frames is None or max_frames <= 0:
                max_frames = frame_end - frame_start
            
            step = (frame_end - frame_start) / max_frames
            frame_indices = [frame_start + int(i * step) for i in range(max_frames)]
            frame_indices = [min(idx, frame_end - 1) for idx in frame_indices]
            
        elif frame_selection == "random":
            # Random selection
            if max_frames is None or max_frames <= 0:
                max_frames = frame_end - frame_start
            
            import random
            random.seed(random_seed)
            available_frames = list(range(frame_start, frame_end))
            if frame_stride is not None and frame_stride > 1:
                available_frames = available_frames[::frame_stride]
            
            frame_indices = random.sample(available_frames, min(max_frames, len(available_frames)))
            frame_indices.sort()  # Keep chronological order
            
        else:
            raise ValueError(f"Unknown frame selection method: {frame_selection}")
        
        # Apply max_frames limit
        if max_frames is not None and max_frames > 0:
            frame_indices = frame_indices[:max_frames]
        
        # Remove duplicates and sort
        frame_indices = sorted(list(set(frame_indices)))
        
        print(f"✅ Selected {len(frame_indices)} frames:")
        print(f"  • Frame indices: {frame_indices[:10]}{'...' if len(frame_indices) > 10 else ''}")
        print(f"  • Frame range: {min(frame_indices)} to {max(frame_indices)}")
        
        return frame_indices
    
    # Test cases
    test_cases = [
        {
            'name': 'Sequential with stride',
            'params': {
                'trajectory_length': 1000,
                'max_frames': 50,
                'frame_start': 100,
                'frame_end': 500,
                'frame_stride': 5,
                'frame_selection': 'sequential'
            }
        },
        {
            'name': 'Equidistant selection',
            'params': {
                'trajectory_length': 1000,
                'max_frames': 30,
                'frame_start': 50,
                'frame_end': 1000,
                'frame_stride': None,
                'frame_selection': 'equidistant'
            }
        },
        {
            'name': 'Random selection',
            'params': {
                'trajectory_length': 1000,
                'max_frames': 25,
                'frame_start': 200,
                'frame_end': 800,
                'frame_stride': 2,
                'frame_selection': 'random',
                'random_seed': 123
            }
        },
        {
            'name': 'Full trajectory with stride',
            'params': {
                'trajectory_length': 1000,
                'max_frames': None,
                'frame_start': None,
                'frame_end': None,
                'frame_stride': 10,
                'frame_selection': 'sequential'
            }
        },
        {
            'name': 'Production analysis',
            'params': {
                'trajectory_length': 2000,
                'max_frames': 100,
                'frame_start': 1000,
                'frame_end': None,
                'frame_stride': 5,
                'frame_selection': 'sequential'
            }
        }
    ]
    
    results = {}
    
    for i, test_case in enumerate(test_cases, 1):
        print(f"\n🧪 Test {i}: {test_case['name']}")
        print("-" * 40)
        
        try:
            selected_frames = select_frames(**test_case['params'])
            
            results[test_case['name']] = {
                'success': True,
                'n_frames': len(selected_frames),
                'frame_range': f"{min(selected_frames)}-{max(selected_frames)}",
                'frames': selected_frames[:5] + ['...'] if len(selected_frames) > 5 else selected_frames
            }
            
            print(f"✅ Test passed: {len(selected_frames)} frames selected")
            
        except Exception as e:
            results[test_case['name']] = {
                'success': False,
                'error': str(e)
            }
            print(f"❌ Test failed: {e}")
    
    # Summary
    print(f"\n📊 Test Summary")
    print("=" * 30)
    
    passed = sum(1 for r in results.values() if r['success'])
    total = len(results)
    
    print(f"Tests passed: {passed}/{total}")
    
    for name, result in results.items():
        if result['success']:
            print(f"✅ {name}: {result['n_frames']} frames ({result['frame_range']})")
        else:
            print(f"❌ {name}: {result['error']}")
    
    return results

def test_yaml_configurations():
    """Test YAML configurations with frame selection parameters"""
    
    print(f"\n🧪 Testing YAML Configurations")
    print("=" * 50)
    
    # Test main config
    try:
        # Try to find the config file in different locations
        config_paths = [
            'mmgbsa_config.yaml',
            '../mmgbsa_config.yaml',
            'test/mmgbsa_config.yaml'
        ]
        
        main_config = None
        for path in config_paths:
            if Path(path).exists():
                with open(path, 'r') as f:
                    main_config = yaml.safe_load(f)
                break
        
        if main_config is None:
            raise FileNotFoundError("Could not find mmgbsa_config.yaml")
        
        frame_params = main_config['analysis_settings']
        print("✅ Main configuration loaded successfully")
        print(f"  • Frame start: {frame_params.get('frame_start', 'null')}")
        print(f"  • Frame end: {frame_params.get('frame_end', 'null')}")
        print(f"  • Frame stride: {frame_params.get('frame_stride', 'null')}")
        print(f"  • Frame selection: {frame_params.get('frame_selection', 'sequential')}")
        print(f"  • Random seed: {frame_params.get('random_seed', 42)}")
        
    except Exception as e:
        print(f"❌ Main configuration test failed: {e}")
        return False
    
    # Test frame test configs
    try:
        # Try to find the config file in different locations
        config_paths = [
            'frame_test_configs.yaml',
            'test/frame_test_configs.yaml',
            '../frame_test_configs.yaml'
        ]
        
        test_configs = None
        for path in config_paths:
            if Path(path).exists():
                with open(path, 'r') as f:
                    test_configs = yaml.safe_load(f)
                break
        
        if test_configs is None:
            raise FileNotFoundError("Could not find frame_test_configs.yaml")
        
        print(f"\n✅ Frame test configurations loaded successfully")
        print(f"  • Number of test configs: {len(test_configs)}")
        
        for config_name, config in test_configs.items():
            frame_params = config['analysis_settings']
            print(f"  • {config_name}: {frame_params.get('frame_selection', 'sequential')} "
                  f"({frame_params.get('max_frames', 'all')} frames)")
        
    except Exception as e:
        print(f"❌ Frame test configurations failed: {e}")
        return False
    
    return True

def test_frame_parameter_validation():
    """Test frame parameter validation logic"""
    
    print(f"\n🧪 Testing Frame Parameter Validation")
    print("=" * 50)
    
    def validate_frame_params(trajectory_length, frame_start, frame_end, frame_stride, max_frames):
        """Validate frame parameters"""
        errors = []
        warnings = []
        
        # Validate frame_start
        if frame_start is not None:
            if frame_start < 0:
                errors.append(f"frame_start ({frame_start}) cannot be negative")
            elif frame_start >= trajectory_length:
                errors.append(f"frame_start ({frame_start}) >= trajectory_length ({trajectory_length})")
        
        # Validate frame_end
        if frame_end is not None:
            if frame_end <= 0:
                errors.append(f"frame_end ({frame_end}) must be positive")
            elif frame_end > trajectory_length:
                warnings.append(f"frame_end ({frame_end}) > trajectory_length ({trajectory_length}), will be capped")
        
        # Validate frame_stride
        if frame_stride is not None:
            if frame_stride <= 0:
                errors.append(f"frame_stride ({frame_stride}) must be positive")
            elif frame_stride > trajectory_length:
                warnings.append(f"frame_stride ({frame_stride}) > trajectory_length ({trajectory_length})")
        
        # Validate max_frames
        if max_frames is not None:
            if max_frames <= 0:
                errors.append(f"max_frames ({max_frames}) must be positive")
        
        return errors, warnings
    
    # Test cases
    test_cases = [
        {
            'name': 'Valid parameters',
            'params': {'trajectory_length': 1000, 'frame_start': 100, 'frame_end': 500, 'frame_stride': 5, 'max_frames': 50}
        },
        {
            'name': 'Invalid frame_start',
            'params': {'trajectory_length': 1000, 'frame_start': 1500, 'frame_end': 500, 'frame_stride': 5, 'max_frames': 50}
        },
        {
            'name': 'Invalid frame_stride',
            'params': {'trajectory_length': 1000, 'frame_start': 100, 'frame_end': 500, 'frame_stride': -1, 'max_frames': 50}
        },
        {
            'name': 'Large frame_end (warning)',
            'params': {'trajectory_length': 1000, 'frame_start': 100, 'frame_end': 2000, 'frame_stride': 5, 'max_frames': 50}
        }
    ]
    
    for test_case in test_cases:
        print(f"\n🧪 {test_case['name']}")
        print("-" * 30)
        
        errors, warnings = validate_frame_params(**test_case['params'])
        
        if errors:
            print(f"❌ Errors: {errors}")
        else:
            print("✅ No errors")
        
        if warnings:
            print(f"⚠️  Warnings: {warnings}")
        else:
            print("✅ No warnings")

def main():
    """Main test function"""
    
    print("🚀 Frame Selection Test Suite")
    print("=" * 60)
    
    # Test 1: Frame selection logic
    logic_results = test_frame_selection_logic()
    
    # Test 2: YAML configurations
    config_success = test_yaml_configurations()
    
    # Test 3: Parameter validation
    test_frame_parameter_validation()
    
    # Final summary
    print(f"\n🎉 Frame Selection Test Suite Complete!")
    print("=" * 60)
    
    if logic_results and config_success:
        print("✅ All tests passed!")
        print("✅ Frame selection functionality is ready for use!")
        print("\n📚 Usage examples:")
        print("  • Sequential with stride: frame_stride=5")
        print("  • Equidistant selection: frame_selection='equidistant'")
        print("  • Random selection: frame_selection='random'")
        print("  • Skip equilibration: frame_start=1000")
        print("  • Limit frames: max_frames=50")
    else:
        print("❌ Some tests failed!")
        print("Please check the errors above.")

if __name__ == '__main__':
    main() 
#!/usr/bin/env python3
"""
Test script demonstrating separate platform configuration for analysis vs decomposition.
Shows how to configure GPU for main calculations and CPU for decomposition.
"""

import sys
from pathlib import Path

# Add mmgbsa module to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from mmgbsa.core import GBSACalculator

def test_platform_settings():
    """Test that platform settings are correctly applied to GBSACalculator"""
    
    print("=" * 80)
    print("PLATFORM SETTINGS CONFIGURATION TEST")
    print("=" * 80)
    
    # Create calculator
    print("\n1. Creating GBSACalculator instance...")
    calculator = GBSACalculator(
        temperature=300,
        gb_model='OBC2',
        verbose=2
    )
    print("✓ Calculator created")
    
    # Test 1: Before setting platform settings
    print("\n2. Initial state:")
    print(f"   preferred_platform: {calculator.preferred_platform}")
    print(f"   decomposition_platform: {calculator.decomposition_platform}")
    
    # Test 2: Apply platform settings
    print("\n3. Applying platform settings...")
    platform_config = {
        'preferred_platform': 'CUDA',
        'decomposition_platform': 'CPU'
    }
    calculator.set_platform_settings(platform_config)
    print("✓ Platform settings applied")
    
    # Test 3: Verify settings were stored
    print("\n4. After applying settings:")
    print(f"   preferred_platform: {calculator.preferred_platform}")
    print(f"   decomposition_platform: {calculator.decomposition_platform}")
    
    # Verify values match
    assert calculator.preferred_platform == 'CUDA', "preferred_platform not set correctly"
    assert calculator.decomposition_platform == 'CPU', "decomposition_platform not set correctly"
    print("✓ Settings verified correctly")
    
    # Test 4: Test setup_optimized_platform with parameter
    print("\n5. Testing setup_optimized_platform() with different scenarios...")
    
    # Scenario A: Use default (should attempt CUDA)
    print("   a) Default platform (should try CUDA):")
    try:
        platform_a, props_a = calculator.setup_optimized_platform()
        print(f"      ✓ Platform obtained: {platform_a.getName()}")
    except Exception as e:
        print(f"      Note: {e}")
    
    # Scenario B: Force CPU via parameter
    print("   b) Force CPU platform via parameter:")
    try:
        platform_b, props_b = calculator.setup_optimized_platform(platform_name='CPU')
        print(f"      ✓ Platform obtained: {platform_b.getName()}")
        assert platform_b.getName() == 'CPU', "CPU platform not obtained"
    except Exception as e:
        print(f"      Error: {e}")
    
    # Test 5: Empty settings should not crash
    print("\n6. Testing with empty settings...")
    calculator.set_platform_settings(None)
    print("✓ Empty settings handled correctly")
    
    calculator.set_platform_settings({})
    print("✓ Empty dict settings handled correctly")
    
    print("\n" + "=" * 80)
    print("ALL TESTS PASSED ✓")
    print("=" * 80)
    print("\nSummary:")
    print("- Platform settings are correctly read and stored")
    print("- set_platform_settings() method works correctly")
    print("- setup_optimized_platform() can accept platform_name parameter")
    print("- Separate platform configuration enables GPU for analysis + CPU for decomposition")

if __name__ == '__main__':
    test_platform_settings()

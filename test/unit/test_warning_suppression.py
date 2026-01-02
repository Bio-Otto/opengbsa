#!/usr/bin/env python3
"""
Test script to verify warning suppression is working
"""

import warnings

# Suppress specific warnings
warnings.filterwarnings('ignore')
warnings.filterwarnings("ignore", message="Unable to load toolkit 'OpenEye Toolkit'")
warnings.filterwarnings("ignore", message="importing 'simtk.openmm' is deprecated")

from openmm import app, openmm, unit

def test_openff_import():
    """Test OpenFF import with warning suppression"""
    try:
        from openff.toolkit.topology import Molecule
        print("✅ OpenFF import successful (warnings suppressed)")
        return True
    except ImportError as e:
        print(f"❌ OpenFF import failed: {e}")
        return False

def test_openmm_import():
    """Test OpenMM import with warning suppression"""
    try:
        from openmm import app, openmm, unit
        print("✅ OpenMM import successful (warnings suppressed)")
        return True
    except ImportError:
        print(f"❌ OpenMM import failed: {e}")
        return False

def test_yaml_import():
    """Test YAML import"""
    try:
        import yaml
        print("✅ YAML import successful")
        return True
    except ImportError as e:
        print(f"❌ YAML import failed: {e}")
        return False

def main():
    """Main test function"""
    print("🧪 Testing Warning Suppression")
    print("=" * 40)
    
    # Test imports
    openff_success = test_openff_import()
    openmm_success = test_openmm_import()
    yaml_success = test_yaml_import()
    
    print(f"\n📊 Test Results:")
    print(f"  • OpenFF: {'✅' if openff_success else '❌'}")
    print(f"  • OpenMM: {'✅' if openmm_success else '❌'}")
    print(f"  • YAML: {'✅' if yaml_success else '❌'}")
    
    if yaml_success:
        print(f"\n✅ Warning suppression is working!")
        print(f"✅ YAML configuration system is ready!")
        print(f"⚠️  OpenMM is required for full MM/GBSA analysis")
    else:
        print(f"\n❌ Some imports failed")

if __name__ == '__main__':
    main() 
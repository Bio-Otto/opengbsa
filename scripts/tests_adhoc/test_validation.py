#!/usr/bin/env python3
"""
Test script for topology validation system.

This script demonstrates the validation system by checking the 6t1h_1151 topology
that has SCNB=1.0 instead of the standard SCNB=2.0.
"""

import sys
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from pathlib import Path
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent))

from mmgbsa.validation import TopologyValidator,
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
 ValidationReport
import parmed
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
 as pmd

def test_6t1h_1151_validation():
    """Test validation on 6t1h_1151 system with SCNB=1.0"""
    
    print("="*70)
    print("TOPOLOGY VALIDATION TEST: 6t1h_1151")
    print("="*70)
    
    # Load topology
    topology_file = 'test/data/6t1h_1151_comp/complex.prmtop'
    
    print(f"\nLoading topology: {topology_file}")
    structure = pmd.load_file(topology_file)
    print(f"✓ Loaded {len(structure.atoms)} atoms")
    
    # Create validation report
    report = ValidationReport()
    
    # Run SCNB/SCEE validation
    print("\n" + "-"*70)
    print("Running SCNB/SCEE validation...")
    print("-"*70)
    scnb_warnings = TopologyValidator.validate_scnb_scee(structure)
    report.add_warnings(scnb_warnings)
    
    # Run force field validation
    print("\n" + "-"*70)
    print("Running force field validation...")
    print("-"*70)
    ff_warnings = TopologyValidator.validate_force_field(structure)
    report.add_warnings(ff_warnings)
    
    # Print summary
    report.print_summary()
    
    # Test result sanity validation
    print("\n" + "="*70)
    print("RESULT SANITY VALIDATION TEST")
    print("="*70)
    
    # Simulate baseline results
    baseline_results = {
        'binding_energy': -28.12,
        'std_dev': 0.43,
        'components': {
            'vdw': -35.90,
            'ele': -3.46,
            'gb': 21.99,
            'sa': -10.75
        }
    }
    
    print(f"\nTesting with baseline results:")
    print(f"  Binding: {baseline_results['binding_energy']:.2f} ± {baseline_results['std_dev']:.2f}")
    print(f"  VDW: {baseline_results['components']['vdw']:.2f}")
    print(f"  ELE: {baseline_results['components']['ele']:.2f}")
    print(f"  GB: {baseline_results['components']['gb']:.2f}")
    print(f"  SA: {baseline_results['components']['sa']:.2f}")
    
    result_report = ValidationReport()
    sanity_warnings = TopologyValidator.validate_system_sanity(
        binding_energy=baseline_results['binding_energy'],
        std_dev=baseline_results['std_dev'],
        components=baseline_results['components']
    )
    result_report.add_warnings(sanity_warnings)
    result_report.print_summary()
    
    # Test with bad results
    print("\n" + "="*70)
    print("TESTING WITH UNUSUAL RESULTS (should trigger warnings)")
    print("="*70)
    
    bad_results = {
        'binding_energy': -133.57,  # From failed SCNB fix
        'std_dev': 25.0,
        'components': {
            'vdw': -141.35,
            'ele': -3.46,
            'gb': 21.99,
            'sa': -10.75
        }
    }
    
    print(f"\nTesting with unusual results:")
    print(f"  Binding: {bad_results['binding_energy']:.2f} ± {bad_results['std_dev']:.2f}")
    
    bad_report = ValidationReport()
    bad_warnings = TopologyValidator.validate_system_sanity(
        binding_energy=bad_results['binding_energy'],
        std_dev=bad_results['std_dev'],
        components=bad_results['components']
    )
    bad_report.add_warnings(bad_warnings)
    bad_report.print_summary()
    
    # Export to dict
    print("\n" + "="*70)
    print("VALIDATION REPORT EXPORT")
    print("="*70)
    
    import json
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

    export = report.to_dict()
    print(json.dumps(export, indent=2))
    
    print("\n" + "="*70)
    print("TEST COMPLETE")
    print("="*70)
    print(f"\nSummary:")
    print(f"  Total warnings: {len(report.warnings)}")
    print(f"  Errors: {sum(1 for w in report.warnings if w.severity == 'ERROR')}")
    print(f"  Warnings: {sum(1 for w in report.warnings if w.severity == 'WARNING')}")
    print(f"  Has errors: {report.has_errors()}")
    print(f"  Has warnings: {report.has_warnings()}")

if __name__ == '__main__':
    test_6t1h_1151_validation()

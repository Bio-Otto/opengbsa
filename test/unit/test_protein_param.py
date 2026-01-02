#!/usr/bin/env python3
"""Direct test of protein parametrization to see full error"""

import sys
sys.path.insert(0, '/home/bio-otto/Desktop/ATHENA-BACKUP/Desktop/mmpbsa/mmgbsa_v0.0.4')

from mmgbsa.core import FixedEnhancedTrueForceFieldMMGBSA

# Create calculator with CHARMM
calculator = FixedEnhancedTrueForceFieldMMGBSA(
    temperature=300,
    verbose=1,
    gb_model='OBC2',
    salt_concentration=0.15,
    use_cache=False,  # Disable cache to see real error
    protein_forcefield='charmm'
)

# Try to parametrize protein
try:
    protein_system, protein_top, protein_pos = calculator.parameterize_protein_amber(
        'analysis_6xj3/6xj3_mdtraj_clean.pdb',
        'UNL'
    )
    print("✓ Success!")
except Exception as e:
    import traceback
    print(f"✗ Error: {e}")
    print("\nFull traceback:")
    traceback.print_exc()

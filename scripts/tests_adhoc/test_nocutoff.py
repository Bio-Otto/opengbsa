#!/usr/bin/env python3
"""
Quick diagnostic: Verify NonbondedForce actually uses NoCutoff
"""
import sys
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

sys.path.insert(0, '/home/bio-otto/Desktop/opengbsa')

import parmed
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
 as pmd
from mmgbsa.mmgbsa_core import StructureManager
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from openmm import app,
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
 unit
import openmm
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


# Load structure
complex_struct = pmd.load_file('test/data/6t1h_1151_comp/complex.prmtop')

# Create system using OpenGBSA method
print("Creating system with implicitSolvent=GBn...")
system = StructureManager.create_openmm_system(
    complex_struct,
    implicitSolvent=app.GBn,
    implicitSolventSaltConc=0.15 * unit.molar
)

# Find NonbondedForce
nb_force = None
for force in system.getForces():
    if isinstance(force, openmm.NonbondedForce):
        nb_force = force
        break

if nb_force:
    method = nb_force.getNonbondedMethod()
    method_names = {0: "NoCutoff", 1: "CutoffNonPeriodic", 2: "CutoffPeriodic", 
                    3: "Ewald", 4: "PME", 5: "LJPME"}
    print(f"✓ NonbondedForce found")
    print(f"  Method: {method} ({method_names.get(method, 'Unknown')})")
    
    if method == 0:
        print("  ✓ CORRECT: Using NoCutoff")
    else:
        print(f"  ❌ PROBLEM: Using cutoffs!")
        print(f"     Cutoff: {nb_force.getCutoffDistance()}")
else:
    print("❌ No NonbondedForce found!")

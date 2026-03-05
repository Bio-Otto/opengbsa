
import parmed
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
 as pmd
import openmm
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from openmm import app,
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
 unit
import sys
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
 as np
sys.path.insert(0, '/home/bio-otto/Desktop/opengbsa')
from mmgbsa.core import StructureManager
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def clean_energy(system, context):
    return context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole)

print("Loading topology & coordinates...")
prmtop_path = 'test/data/6t1h_1151_comp/complex.prmtop'
pdb_path = 'test/data/6t1h_1151_comp/complex.pdb'

complex_struct = pmd.load_file(prmtop_path)
pdb = pmd.load_file(pdb_path)

if len(complex_struct.atoms) != len(pdb.atoms):
    print("Error: Atom count mismatch between prmtop and pdb")
    sys.exit(1)

complex_struct.coordinates = pdb.coordinates
complex_struct.box = pdb.box

cutoffs = [12.0, 999.0]
print(f"{'Cutoff (A)':<15} {'VDW+Elec (kcal/mol)':<25}")
print("-" * 40)

for cut in cutoffs:
    if cut >= 100.0:
        method = app.NoCutoff
        cut_str = "NoCutoff"
        dist = 999.0*unit.angstrom
    else:
        method = app.CutoffNonPeriodic
        cut_str = f"{cut}"
        dist = cut*unit.angstrom
        
    try:
        # Create system WITHOUT implicit solvent to avoid GBSA forces complicating things
        # But we still want to test the force split if it happens with implicitSolvent=None?
        # No, creating with implicitSolvent=None keeps NonbondedForce intact usually.
        # But we want to test exactly what OpenGBSA does (with ImplicitSolvent=GBn).
        
        system = StructureManager.create_openmm_system(
            complex_struct,
            implicitSolvent=app.GBn, 
            implicitSolventSaltConc=0.15*unit.molar,
            nonbondedMethod=method
        )
        
        # Remove unwanted forces to isolate VDW/Elec and avoid cutoff conflicts
        # Keep NonbondedForce (Group 0) and CustomNonbondedForce (Group 3)
        # Remove GBSAOBCForce, CustomGBForce, CMMotionRemover, etc.
        
        forces_to_keep = []
        for f in system.getForces():
            if isinstance(f, (openmm.NonbondedForce, openmm.CustomNonbondedForce)):
                forces_to_keep.append(f)
                
                # Enforce cutoff
                if hasattr(f, 'setCutoffDistance'):
                    f.setCutoffDistance(dist)
                if hasattr(f, 'setNonbondedMethod'):
                     try:
                        if isinstance(f, openmm.NonbondedForce):
                            f.setNonbondedMethod(method)
                        elif isinstance(f, openmm.CustomNonbondedForce):
                            # CustomNonbondedForce enums might differ
                            if method == app.NoCutoff: f.setNonbondedMethod(openmm.CustomNonbondedForce.NoCutoff)
                            elif method == app.CutoffNonPeriodic: f.setNonbondedMethod(openmm.CustomNonbondedForce.CutoffNonPeriodic)
                     except:
                        pass
        
        # Clear forces and re-add kept ones
        # Use a new system or modify existing? simpler to modify
        # But getForces returns copy? No, reference.
        # removeForce by index. Iterate backwards.
        
        for i in range(system.getNumForces())[::-1]:
            f = system.getForce(i)
            if not isinstance(f, (openmm.NonbondedForce, openmm.CustomNonbondedForce)):
                system.removeForce(i)
                
        # Calculate
        integrator = openmm.VerletIntegrator(0.001)
        context = openmm.Context(system, integrator)
        context.setPositions(complex_struct.coordinates) 
        
        e_pot = clean_energy(system, context)
        
        print(f"{cut_str:<15} {e_pot:<25.2f}")
        
    except Exception as e:
        print(f"{cut_str:<15} ERROR: {e}")

print("-" * 40)


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
from mmgbsa.core import StructureManager,
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
 GBSACalculator

def test_lcpo_energy():
    print("Testing LCPO Surface Area Energy...")
    
    # Load inputs
    complex_prmtop = '/home/bio-otto/Desktop/opengbsa/test/data/6t1h_1151_comp/complex.prmtop'
    complex_pdb = '/home/bio-otto/Desktop/opengbsa/test/data/6t1h_1151_comp/complex.pdb'
    
    struct = pmd.load_file(complex_prmtop, xyz=complex_pdb)
    
    # Initialize Manager with LCPO
    manager = StructureManager()
    manager.load_structure(struct)
    manager.sa_model = 'LCPO'  # Force LCPO model
    
    # Create System
    # Note: implicitSolvent='OBC2' will trigger _setup_implicit_solvent which calls _setup_surface_area_force
    try:
        system = manager.create_openmm_system(
            implicitSolvent='OBC2',
            nonbondedMethod='NoCutoff',
            solute_dielectric=1.0,
            solvent_dielectric=78.5,
            salt_concentration=0.15 * unit.molar
        )
    except Exception as e:
        print(f"Error creating system: {e}")
        import traceback
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

        traceback.print_exc()
        return

    # Create Context
    integrator = openmm.VerletIntegrator(1.0 * unit.femtoseconds)
    context = openmm.Context(system, integrator)
    context.setPositions(struct.positions)
    
    # Calculate Energy
    state = context.getState(getEnergy=True)
    total_energy = state.getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole)
    print(f"Total Potential Energy: {total_energy:.4f} kcal/mol")
    
    # Enhance: Decompose to find SA term
    # We iterate forces to find LCPOForce
    sa_energy = 0.0
    for force in system.getForces():
        name = force.__class__.__name__
        if 'LCPOForce' in name:
            # Create a temp context for just this force to be sure, 
            # or just rely on group if we knew it separately.
            # But simpler: use Core's decompose method/logic logic manually here
            # Actually, let's just isolate this force
            
            # Temporary system for SA
            sa_system = openmm.System()
            for i in range(system.getNumParticles()):
                sa_system.addParticle(system.getParticleMass(i))
            
            # Copy force (shallow copy might work, but safer to add the same object reference if openmm allows)
            # OpenMM forces are objects, can be in multiple systems? No, usually not.
            # But we can read the energy from the main context if we used force groups.
            # In core.py, _setup_lcpo_force sets group 4.
            
            group = force.getForceGroup()
            print(f"Found {name} in Group {group}")
            
            state_sa = context.getState(getEnergy=True, groups={group})
            sa_energy = state_sa.getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole)
            print(f"Surface Area Energy (LCPO): {sa_energy:.4f} kcal/mol")
            break
    else:
        print("Error: LCPOForce not found in system!")

if __name__ == "__main__":
    test_lcpo_energy()

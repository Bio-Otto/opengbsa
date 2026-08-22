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

sys.path.insert(0, '/home/bio-otto/Desktop/opengbsa')
from mmgbsa.mmgbsa_core import StructureManager
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


print("Loading complex.prmtop...")
struct = pmd.load_file('test/data/6t1h_1151_comp/complex.prmtop')

print("Creating System with openmm.app.GBn...")
try:
    system = StructureManager.create_openmm_system(
        struct, 
        implicitSolvent=None,
        implicitSolventSaltConc=0.0*unit.molar,
        nonbondedMethod=app.CutoffNonPeriodic
    )
    
    # Check NonbondedForce
    forces = [f for f in system.getForces() if isinstance(f, openmm.NonbondedForce)]
    if forces:
        nb = forces[0]
        print(f"NonbondedForce found with {nb.getNumParticles()} particles")
        c, s, e = nb.getParticleParameters(0)
        print(f"Atom 0 (NB): Eps={e}")
    
    # Check CustomNonbondedForce
    custom_forces = [f for f in system.getForces() if isinstance(f, openmm.CustomNonbondedForce)]
    if custom_forces:
        print(f"CustomNonbondedForce found!")
        cf = custom_forces[0]
        print(f"Energy function: {cf.getEnergyFunction()}")
        if cf.getNumParticles() > 0:
            params = cf.getParticleParameters(0)
            print(f"Atom 0 (Custom) params: {params}")

        
except Exception as e:
    print(f"Error: {e}")


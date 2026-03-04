
import sys
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import os
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from unittest.mock import MagicMock
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import types
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


# Create a dummy module for plotly
plotly = types.ModuleType("plotly")
plotly.graph_objects = MagicMock()
plotly.express = MagicMock()
plotly.io = MagicMock()

sys.modules["plotly"] = plotly
sys.modules["plotly.graph_objects"] = plotly.graph_objects
sys.modules["plotly.express"] = plotly.express
sys.modules["plotly.io"] = plotly.io

# Mock openff
openff = types.ModuleType("openff")
openff.toolkit = types.ModuleType("openff.toolkit")
openff.toolkit.topology = types.ModuleType("openff.toolkit.topology")
openff.toolkit.typing = types.ModuleType("openff.toolkit.typing")
openff.toolkit.typing.engines = types.ModuleType("openff.toolkit.typing.engines")
openff.toolkit.typing.engines.smirnoff = types.ModuleType("openff.toolkit.typing.engines.smirnoff")

# Mock classes
class MockMolecule:
    pass
class MockForceField:
    pass

openff.toolkit.topology.Molecule = MockMolecule
openff.toolkit.typing.engines.smirnoff.ForceField = MockForceField

sys.modules["openff"] = openff
sys.modules["openff.toolkit"] = openff.toolkit
sys.modules["openff.toolkit.topology"] = openff.toolkit.topology
sys.modules["openff.toolkit.typing"] = openff.toolkit.typing
sys.modules["openff.toolkit.typing.engines"] = openff.toolkit.typing.engines
sys.modules["openff.toolkit.typing.engines.smirnoff"] = openff.toolkit.typing.engines.smirnoff

# Mock openmmforcefields
openmmforcefields = types.ModuleType("openmmforcefields")
openmmforcefields.generators = types.ModuleType("openmmforcefields.generators")

class MockSystemGenerator:
    pass

openmmforcefields.generators.SystemGenerator = MockSystemGenerator

sys.modules["openmmforcefields"] = openmmforcefields
sys.modules["openmmforcefields.generators"] = openmmforcefields.generators

sys.path.append(os.getcwd())

sys.path.append(os.getcwd())
from mmgbsa.core import GBSACalculator,
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
 StructureManager
from openmm import app,
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
 unit, openmm
import parmed
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
 as pmd
import numpy
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
 as np

def check_system_params():
    print("--- Debugging System Parameters ---")
    
    # Paths
    complex_prmtop = "test/data/6t1h_1151_comp/complex.prmtop"
    complex_pdb_path = "test/data/6t1h_1151_comp/temp_fixed_for_panda.pdb" # Use PDB for coords if needed
    
    # 1. Load Prmtop directly with ParmEd
    print(f"Loading {complex_prmtop} with ParmEd...")
    struct = pmd.load_file(complex_prmtop)
    
    # Extract Radii from Prmtop
    # Amber prmtop stores radii in .radii usually, or we can get via .sigma calculation?
    # GB radii are usually in 'RADII' section if mbondi2 is used.
    # ParmEd stores them in atom.radii if available, or atom.solvent_radius?
    # Let's check a few atoms.
    print(f"Prmtop Atoms: {len(struct.atoms)}")
    prmtop_radii = [a.solvent_radius for a in struct.atoms]
    print(f"Sample Prmtop Radii (Res 0-2): {prmtop_radii[:10]}")
    
    # 2. Create System via GBSACalculator logic
    # We need to simulate how core.py creates it.
    calc = GBSACalculator(gb_model='GBn', sa_model='LCPO', salt_concentration=0.15)
    
    # Mock systems_dict logic from core.py (prepare_decomposition_systems or build_complex_system)
    # core.py: build_complex_system(complex_pdb, ...)
    # It uses StructureManager.load_complex -> create_openmm_system
    
    print("\nCreating OpenMM System via GBSACalculator logic...")
    # We replicate the steps from core.py roughly
    
    # Step A: Load Structure
    # If we pass prmtop to load_complex, it loads it.
    complex_structure = StructureManager.load_complex(complex_pdb_path, prmtop_path=complex_prmtop)
    
    # Step B: Create System
    # core.py: create_openmm_system(structure, implicitSolvent=app.GBn, ...)
    # Note: refined_gbsa_forces is called AFTER.
    
    system = complex_structure.createSystem(
        nonbondedMethod=app.NoCutoff,
        constraints=None,
        implicitSolvent=app.GBn,
        implicitSolventSaltConc=0.15*unit.molar
    )
    
    # Step C: Refine Forces (This is where Radii might get messed up)
    calc.gbsa_manager.refine_gbsa_forces(system, complex_structure.topology)
    
    # 3. Inspect OpenMM System Radii
    print("\nInspecting OpenMM System GB Forces...")
    gb_forces = [f for f in system.getForces() if isinstance(f, (openmm.GBSAOBCForce, openmm.CustomGBForce))]
    
    for f in gb_forces:
        print(f"Found GB Force: {type(f).__name__}")
        if isinstance(f, openmm.GBSAOBCForce):
             # Check radii
             # particle parameters: charge, radius, scalingFactor
             n = f.getNumParticles()
             omm_radii = []
             for i in range(min(10, n)):
                 q, r, s = f.getParticleParameters(i)
                 omm_radii.append(r.value_in_unit(unit.angstroms))
             print(f"  OpenMM Radii (First 10, Angstroms): {omm_radii}")
             
             # Compare
             print(f"  Prmtop Radii (First 10, Angstroms): {prmtop_radii[:10]}")
             
             diffs = [abs(o - p) for o, p in zip(omm_radii, prmtop_radii[:10])]
             print(f"  Differences: {diffs}")
             
        elif isinstance(f, openmm.CustomGBForce):
            print("  (CustomGBForce parameters are complex to read directly, skipping detail check for now)")
            
    # 4. Inspect SA Force (LCPO)
    print("\nInspecting SA Force...")
    sa_forces = [f for f in system.getForces() if isinstance(f, openmm.LCPOForce) or (isinstance(f, openmm.GBSAOBCForce) and f.getSurfaceAreaEnergy().value_in_unit(unit.kilojoule_per_mole/unit.nanometer**2) > 0)]
    
    if not sa_forces:
        # Check CustomGB with SA?
        pass
        
    for f in sa_forces:
        print(f"Found SA Force: {type(f).__name__}")
        if isinstance(f, openmm.LCPOForce):
             # Check parameters 
             # addParticle(radius, p1, p2, p3, p4)
             n = f.getNumParticles()
             lcpo_radii = []
             for i in range(min(10, n)):
                 pars = f.getParticleParameters(i)
                 # radius is pars[0]
                 lcpo_radii.append(pars[0].value_in_unit(unit.angstroms))
             print(f"  LCPO Radii (First 10, Angstroms): {lcpo_radii}")
             print(f"  Prmtop Radii (First 10, Angstroms): {prmtop_radii[:10]}")

if __name__ == "__main__":
    check_system_params()

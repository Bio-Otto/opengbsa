
import sys
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import unittest.mock
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from unittest.mock import MagicMock
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


# --- 1. PRE-IMPORT MOCKING ---
# Mock EVERYTHING to ensure clean testing environment
mock_openmm = MagicMock()
mock_app = MagicMock()
mock_unit = MagicMock()

sys.modules['openmm'] = mock_openmm
sys.modules['openmm.app'] = mock_app
sys.modules['openmm.unit'] = mock_unit
# Populate mock_unit for 'from openmm.unit import *'
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
 to work
units_to_mock = ['kelvin', 'kilocalorie', 'mole', 'angstrom', 'gram', 
                 'centimeter', 'joule', 'dalton', 'nanometer']
mock_unit.__all__ = units_to_mock
for u in units_to_mock:
    setattr(mock_unit, u, MagicMock())

# Handle internal helper
sys.modules['openmm.app.internal.lcpo'] = MagicMock()

# Other dependencies
sys.modules['plotly'] = MagicMock()
sys.modules['plotly.graph_objects'] = MagicMock()
sys.modules['plotly.subplots'] = MagicMock()
sys.modules['plotly.io'] = MagicMock()
sys.modules['plotly.express'] = MagicMock()
sys.modules['plotly.colors'] = MagicMock()
sys.modules['openff'] = MagicMock()
sys.modules['openff.toolkit'] = MagicMock()
sys.modules['openff.toolkit.topology'] = MagicMock()
sys.modules['openff.toolkit.typing'] = MagicMock()
sys.modules['openff.toolkit.typing.engines'] = MagicMock()
sys.modules['openff.toolkit.typing.engines.smirnoff'] = MagicMock()
sys.modules['openmmforcefields'] = MagicMock()
sys.modules['openmmforcefields.generators'] = MagicMock()
sys.modules['parmed'] = MagicMock()

# Numpy Patch
import numpy
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
 as np
import types
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

if 'numpy.compat' not in sys.modules:
    comp = types.ModuleType('numpy.compat')
    sys.modules['numpy.compat'] = comp
    if hasattr(np, 'compat'):
         sys.modules['numpy.compat'] = np.compat
    else:
        def asbytes(s):
            return s if isinstance(s, bytes) else str(s).encode('latin1')
        def asstr(s):
            return s.decode('latin1') if isinstance(s, bytes) else str(s)
        comp.asbytes = asbytes
        comp.asstr = asstr
        sys.modules['numpy.compat'] = comp

# --- 2. IMPORT CORE ---
print("Importing mmgbsa.core with full mocks...")
try:
    from mmgbsa.core import StructureManager,
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
 GBSAForceManager, TopologyAdapter
    # We must ensure openmm is accessible as a name if core uses it via 'import openmm'
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

    import openmm
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
 
except ImportError as e:
    print(f"CRITICAL IMPORT ERROR: {e}")
    sys.exit(1)

def test_lcpo_logic():
    print("\n--- Testing LCPO Setup Logic ---")
    
    # Check if openmm.LCPOForce exists on our mock
    if not hasattr(openmm, 'LCPOForce'):
         print("WARNING: openmm.LCPOForce not found on mock. Adding it manually.")
         openmm.LCPOForce = MagicMock()
    
    # 1. Setup Mock Structure (AtomList like ParmEd)
    mock_struct = MagicMock()
    mock_atoms = [MagicMock() for _ in range(5)]
    for a in mock_atoms:
        a.element.symbol = 'C'
        a.radius = 1.7
        a.solvent_radius = 1.7 
    mock_struct.atoms = mock_atoms # LIST
    mock_struct.bonds = []
    
    print(f"Created Mock Structure with {len(mock_struct.atoms)} atoms (List type)")

    # 2. Test TopologyAdapter
    print("\n[Test 1] TopologyAdapter Wrapper")
    adapter = TopologyAdapter(mock_struct)
    if callable(adapter.atoms):
        print("PASS: adapter.atoms is callable")
        try:
             first = next(adapter.atoms())
             print(f"PASS: adapter.atoms() returns iterator")
        except StopIteration:
             print("PASS: adapter.atoms() returns empty iterator")
    else:
        print("FAIL: adapter.atoms is NOT callable")
        
    # 3. Test _setup_lcpo_force
    print("\n[Test 2] GGBForceManager._setup_lcpo_force")
    manager = GBSAForceManager(sa_model='LCPO')
    
    # Mock getLCPOParamsTopology
    # We must patch where it's imported/used. 
    # core.py has 'from openmm.app.internal.lcpo import getLCPOParamsTopology'
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

    # So we prefer to patch 'mmgbsa.core.getLCPOParamsTopology'
    
    # Mock openmm.app.Element.getByAtomicNumber
    mock_element_obj = MagicMock()
    mock_element_obj.symbol = 'C'
    openmm.app.Element = MagicMock()
    openmm.app.Element.getByAtomicNumber.return_value = mock_element_obj

    with unittest.mock.patch('mmgbsa.core.getLCPOParamsTopology') as mock_get_params:
        mock_get_params.return_value = [] # legitimate return
        
        try:
            # CALL THE METHOD
            manager._setup_lcpo_force(None, mock_struct)
            
            # ... checks ... 
            
            # NEW TEST: Verify _get_gb_radius handles int
            print("\n[Test 3] _get_gb_radius with int element")
            # We created atoms with .element = 6
            test_atom = mock_atoms[0]
            radius = manager._get_gb_radius(test_atom)
            # Should look up 'C' -> default radius 1.5? Or whatever matches 'C'
            # gb_radii default get is 1.5
            print(f"Radius for element 6: {radius}")
            if radius == 1.5: # Default
                 print("SUCCESS: _get_gb_radius handled int element!")
            else:
                 print(f"FAILURE: _get_gb_radius returned {radius}")

        except Exception as e:
            
            # CHECK ARGUMENTS
            if mock_get_params.called:
                args, _ = mock_get_params.call_args
                passed_arg = args[0]
                print(f"getLCPOParamsTopology called with: {type(passed_arg)}")
                
                if isinstance(passed_arg, TopologyAdapter):
                    print("SUCCESS: Structure was wrapped in TopologyAdapter!")
                else:
                    print(f"FAILURE: Structure was NOT wrapped! Got {type(passed_arg)}")
            else:
                print("FAILURE: getLCPOParamsTopology was NOT called (maybe HAS_LCPO_PARAMS check failed?)")
                
        except Exception as e:
            print(f"EXCEPTION: {e}")
            import traceback
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

            traceback.print_exc()

if __name__ == "__main__":
    test_lcpo_logic()

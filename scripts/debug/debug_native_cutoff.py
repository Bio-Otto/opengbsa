
import sys
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import os
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from pathlib import Path
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import openmm.app
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
 as app
import openmm
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
 as mm
import openmm.unit
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
 as unit

# Add project root to path
sys.path.insert(0, os.path.abspath('.'))
sys.path.append(os.path.abspath('mmgbsa'))

from mmgbsa.topology import TopologyLoader
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from mmgbsa.inputs import EngineMode
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def inspect_system(cutoff_val, desc):
    print(f"\n--- Inspeting System: {desc} (Cutoff: {cutoff_val}) ---")
    sys.stdout.flush()
    
    complex_prmtop = "test/data/6t1h_6466_comp/complex.prmtop"
    if not Path(complex_prmtop).exists():
        print(f"Error: {complex_prmtop} not found.")
        sys.stdout.flush()
        return

    # Mocking what GBSACalculator.parameterize_protein_amber does in Native Mode
    loader_kwargs = {
        'implicitSolvent': app.OBC2, # Mocking
        'implicitSolventSaltConc': 0.15 * unit.molar,
        'nonbondedMethod': app.NoCutoff
    }
    
    if cutoff_val is not None:
        loader_kwargs['nonbondedMethod'] = app.CutoffNonPeriodic
        loader_kwargs['nonbondedCutoff'] = cutoff_val * unit.angstroms
        print(f"Requesting CutoffNonPeriodic with {cutoff_val} A")
        sys.stdout.flush()
    else:
        print("Requesting NoCutoff")
        sys.stdout.flush()

    try:
        # Load system using correctly detected mode (AMBER)
        system, topology, _ = TopologyLoader.load_system(complex_prmtop, EngineMode.AMBER, **loader_kwargs)
        
        # Inspect Forces
        nb_force = None
        custom_nb_force = None
        
        for f in system.getForces():
            if isinstance(f, mm.NonbondedForce):
                nb_force = f
            elif isinstance(f, mm.CustomNonbondedForce):
                custom_nb_force = f
                
        if nb_force:
            method_idx = nb_force.getNonbondedMethod()
            methods = {0: 'NoCutoff', 1: 'CutoffNonPeriodic', 2: 'CutoffPeriodic', 3: 'Ewald', 4: 'PME'}
            method_name = methods.get(method_idx, str(method_idx))
            
            cutoff_dist = nb_force.getCutoffDistance()
            print(f"NonbondedForce Method: {method_name}")
            print(f"NonbondedForce Cutoff: {cutoff_dist}")
            sys.stdout.flush()
            
            if cutoff_val is not None:
                if method_idx != 1: # CutoffNonPeriodic
                    print("FAILURE: Expected CutoffNonPeriodic (1), got", method_idx)
                if abs(cutoff_dist.value_in_unit(unit.angstroms) - cutoff_val) > 0.001:
                    print(f"FAILURE: Cutoff mismatch! Expected {cutoff_val}, got {cutoff_dist.value_in_unit(unit.angstroms)}")
            else:
                if method_idx != 0: # NoCutoff
                     print("FAILURE: Expected NoCutoff (0), got", method_idx)
            sys.stdout.flush()
                     
        else:
            print("FAILURE: No NonbondedForce found!")
            sys.stdout.flush()

    except Exception as e:
        print(f"Exception during loading: {e}")
        sys.stdout.flush()

if __name__ == "__main__":
    inspect_system(10.0, "Finite Cutoff 10.0")
    inspect_system(12.0, "Finite Cutoff 12.0")
    inspect_system(None, "No Cutoff (Null)")

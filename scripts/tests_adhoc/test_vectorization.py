import numpy
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
 as np
import warnings
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


# Mock classes/functions if needed, or import from
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
 decomposition.py
# But decomposition.py imports openmm, which might fail if env issue.
# I'll just copy the vectorized logic here to test SYNTAX and SHAPE match.
# Or better: try to import the
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
 function.

try:
    from mmgbsa.decomposition import _calculate_pairwise_interactions_standalone
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

    print("Import successful")
except ImportError:
    print("Import failed (expected if dependencies missing in base python)")
    exit(0)

# Mock Data
n_atoms = 100
pos_nm = np.random.rand(n_atoms, 3).astype(np.float32)
residue_map = {'RES_1': [0, 1, 2], 'RES_2': [3, 4, 5]}
ligand_indices = list(range(n_atoms)) # All atoms
exclusions = {(0, 1): (0.0, 0.0, 0.0)} # Bonded exclusion between 0 and 1

# Mock System context? function needs params but we can mock pdb_params
pdb_params = {i: (0.1, 0.3, 0.1) for i in range(n_atoms)} 
# q=0.1, sig=0.3 nm, eps=0.1 kJ/mol

print("Running standalone calc...")
try:
    res = _calculate_pairwise_interactions_standalone(
        system=None, # unused if pdb_params provided?
        context=None, 
        positions=pos_nm, 
        residue_map=residue_map, 
        ligand_indices=ligand_indices, 
        salt_concentration=0.15, 
        pdb_params=pdb_params, 
        exclusions=exclusions
    )
    print("Success!")
    print(res)
except Exception as e:
    import traceback
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

    traceback.print_exc()
    print(f"FAILED: {e}")

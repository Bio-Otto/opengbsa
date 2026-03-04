import os
import sys
sys.path.insert(0, os.getcwd())
from pandamap.core import HybridProtLigMapper
from pandamap.create_3d_view import create_pandamap_3d_viz

output_file = "test_pandamap.html"
# We need a PDB complex. Let's find one.
import glob
pdbs = glob.glob("test/data/*.pdb") + glob.glob("*.pdb") + ["test_dummy.pdb"]
if len(pdbs) > 0 and os.path.exists(pdbs[0]):
    print(f"Using {pdbs[0]}")
    try:
        mapper = HybridProtLigMapper(pdbs[0])
        mapper.calculate_interactions()
        create_pandamap_3d_viz(mapper, output_file)
        print("Done")
    except Exception as e:
        print(f"Error: {e}")
else:
    print("No PDB found")

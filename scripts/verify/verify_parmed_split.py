
import sys
import types
import numpy as np

# Patch numpy.compat for ParmEd (NumPy 2.x compatibility)
if 'numpy.compat' not in sys.modules:
    comp = types.ModuleType('numpy.compat')
    sys.modules['numpy.compat'] = comp
    if hasattr(np, 'compat'):
         sys.modules['numpy.compat'] = np.compat
    else:
        def asbytes(s):
            if isinstance(s, bytes):
                return s
            return str(s).encode('latin1')
            
        def asstr(s):
            if isinstance(s, bytes):
                return s.decode('latin1')
            return str(s)
            
        comp.asbytes = asbytes
        comp.asstr = asstr
        sys.modules['numpy.compat'] = comp

import parmed as pmd
from openmm import app

def test_parmed_split():
    print("--- Testing ParmEd Split Stability ---")
    
    prmtop_path = "test/data/6t1h_1151_comp/complex.prmtop"
    pdb_path = "test/data/6t1h_1151_comp/temp_fixed_for_panda.pdb"
    
    print(f"Loading {prmtop_path}...")
    struct = pmd.load_file(prmtop_path)
    print(f"Original Atoms: {len(struct.atoms)}")
    print(f"Class: {type(struct)}")
    
    # Simulate Step 1: Strip Solvent
    print("\n[Step 1] Stripping Solvent...")
    # mask = ":WAT,HOH,TIP3,SOL,Na+,Cl-,Mg+,K+,Zn+,Ca+,NA,CL,MG,K,ZN,CA"
    # Simplified mask for test
    mask = ":WAT,Na+,Cl-" 
    try:
        struct.strip(mask)
        print(f"Stripped Atoms: {len(struct.atoms)}")
    except Exception as e:
        print(f"Step 1 Failed: {e}")
        return

    # Simulate Step 2: Split Components (The crashing part)
    print("\n[Step 2] Splitting Components (Copy + Strip)...")
    ligand_resname = "LIG"
    
    try:
        target_cls = type(struct)
        print(f"Target Class: {target_cls}")
        
        # Copy for Receptor
        print("Copying for Receptor...")
        receptor = struct.copy(cls=target_cls)
        print(f"Receptor Copy Atoms: {len(receptor.atoms)}")
        
        print(f"Stripping Ligand :{ligand_resname} from Receptor...")
        receptor.strip(f":{ligand_resname}")
        print(f"Receptor Atoms: {len(receptor.atoms)}")
        
        # Copy for Ligand
        print("Copying for Ligand...")
        ligand = struct.copy(cls=target_cls)
        print(f"Stripping Receptor !:{ligand_resname} from Ligand...")
        ligand.strip(f"!:{ligand_resname}")
        print(f"Ligand Atoms: {len(ligand.atoms)}")
        
    except Exception as e:
        print(f"Step 2 Failed (Copy+Strip): {e}")
        import traceback
        traceback.print_exc()

    # Simulate Step 3: Create System (to check GB)
    print("\n[Step 3] Creating OpenMM System (Check GB)...")
    try:
        # Need to verify if createSystem works and produces GB
        # We need a dummy method
        # Implicit solvent needs to be passed?
        # For AmberParm, it uses pointers?
        s = receptor.createSystem(implicitSolvent=app.OBC2)
        print("System created successfully.")
        # Check forces
        forces = [f.__class__.__name__ for f in s.getForces()]
        print(f"Forces: {forces}")
    except Exception as e:
        print(f"Step 3 Failed (CreateSystem): {e}")

    # TEST ALTERNATIVE: Slicing
    print("\n--- Testing SLICING Alternative ---")
    struct2 = pmd.load_file(prmtop_path)
    print("Loaded fresh struct.")
    
    try:
        # Global Strip via Slice
        # Keep NOT (Water/Ions)
        print("Slicing to remove water...")
        dry_struct = struct2['!(:WAT,Na+,Cl-)']
        print(f"Dry Atoms: {len(dry_struct.atoms)}")
        print(f"Dry Class: {type(dry_struct)}")
        
        # Receptor Slice
        print("Slicing Receptor...")
        rec_slice = dry_struct[f'!:{ligand_resname}']
        print(f"Rec Slice Atoms: {len(rec_slice.atoms)}")
        
        # Ligand Slice
        print("Slicing Ligand...")
        lig_slice = dry_struct[f':{ligand_resname}']
        print(f"Lig Slice Atoms: {len(lig_slice.atoms)}")
        
        # Check System Creation on Slice
        print("Creating System from Slice...")
        sys_slice = rec_slice.createSystem(implicitSolvent=app.OBC2)
        print("Slice System created.")
        forces = [f.__class__.__name__ for f in sys_slice.getForces()]
        print(f"Slice Forces: {forces}")
        
    except Exception as e:
        print(f"Slicing Failed: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    test_parmed_split()

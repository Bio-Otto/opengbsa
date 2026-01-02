import sys
import numpy
import types

# Patch numpy.compat for ParmEd
try:
    from numpy import compat
except ImportError:
    comp = types.ModuleType('numpy.compat')
    def asbytes(s): return s if isinstance(s, bytes) else str(s).encode('latin1')
    def asstr(s): return s if isinstance(s, str) else s.decode('latin1')
    comp.asbytes = asbytes
    comp.asstr = asstr
    numpy.compat = comp
    sys.modules['numpy.compat'] = comp

try:
    import parmed
    print(f"ParmEd version: {parmed.__version__}")
    
    # Try explicit format import
    print("Attempting to read TPR with GromacsTopologyFile...")
    try:
        from parmed.gromacs import GromacsTopologyFile
        top = GromacsTopologyFile('test/complex2.tpr')
        print(f"Loaded successfully! Atoms: {len(top.atoms)}")
        
        top.save('test/complex2_ref.pdb', overwrite=True)
        print("Saved test/complex2_ref.pdb")
    except Exception as e:
        print(f"Explicit load failed: {e}")
        
except Exception as e:
    print(f"Import failed: {e}")

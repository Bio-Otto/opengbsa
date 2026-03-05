
import openmm
from openmm import app, unit
import sys

def inspect_gbn():
    prmtop = app.AmberPrmtopFile('test/data/6t1h_1151_comp/complex.prmtop')
    system = prmtop.createSystem(implicitSolvent=app.GBn)
    
    print("Forces in System:")
    for f in system.getForces():
        print(f"  - {f.__class__.__name__}")
        if isinstance(f, openmm.CustomGBForce):
            # CustomGBForce has multiple terms, not a single function getter
            # print(f"    Energy Function: {f.getEnergyFunction()}")
            print("    Global Parameters:")
            for i in range(f.getNumGlobalParameters()):
                name = f.getGlobalParameterName(i)
                val = f.getGlobalParameterDefaultValue(i)
                print(f"      - {name}: {val}")
            
            print("    Energy Terms:")
            for i in range(f.getNumEnergyTerms()):
                expr, type = f.getEnergyTermParameters(i)
                print(f"      - Term {i}: {expr} (Type: {type})")
            
            print("    Per-Particle Parameters:")
            for i in range(f.getNumPerParticleParameters()):
                name = f.getPerParticleParameterName(i)
                print(f"      - {name}")

if __name__ == "__main__":
    inspect_gbn()

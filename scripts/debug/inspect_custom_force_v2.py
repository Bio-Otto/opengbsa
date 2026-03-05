
import openmm
import openmm.app as app
import sys

# Redirect output to file to ensure we capture it even if buffered
sys.stdout = open('inspect_force_v2.log', 'w')

print("Loading prmtop...")
try:
    prmtop = app.AmberPrmtopFile('test/data/6t1h_1151_comp/complex.prmtop')
    print("Creating system...")
    # Matches core.py usage: nonbondedMethod=NoCutoff, implicitSolvent=OBC2 (approx)
    # core.py uses: implicitSolvent=app.OBC2 (which is string 'OBC2' mapped to method)
    # Let's check what createSystem produces by default for Amber
    system = prmtop.createSystem(nonbondedMethod=app.NoCutoff, implicitSolvent=app.OBC2)
    
    print("Forces in system:")
    for i, f in enumerate(system.getForces()):
        print(f"Force {i}: {f.__class__.__name__}")
        if isinstance(f, openmm.CustomNonbondedForce):
            print(f"  Energy Expression: {f.getEnergyFunction()}")
            print(f"  Num Particles: {f.getNumParticles()}")
            print(f"  Num Global Params: {f.getNumGlobalParameters()}")
            for k in range(f.getNumGlobalParameters()):
                 print(f"    Param {k}: {f.getGlobalParameterName(k)} = {f.getGlobalParameterDefaultValue(k)}")
            
            # Check interaction groups?
            print(f"  Num Interaction Groups: {f.getNumInteractionGroups()}")
            
            # Check first few particles
            if f.getNumParticles() > 0:
                 print(f"    Particle 0 params: {f.getParticleParameters(0)}")

        if isinstance(f, openmm.GBSAOBCForce):
            print(f"  GBSAOBCForce found. SoluteDielectric: {f.getSoluteDielectric()}, SolventDielectric: {f.getSolventDielectric()}")
            
except Exception as e:
    print(f"Error: {e}")


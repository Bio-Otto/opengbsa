
import openmm
import openmm.app as app
import sys

# Load prmtop
prmtop = app.AmberPrmtopFile('test/data/6t1h_1151_comp/complex.prmtop')
system = prmtop.createSystem(nonbondedMethod=app.NoCutoff)

print("Forces in system:")
for i, f in enumerate(system.getForces()):
    print(f"Force {i}: {f.__class__.__name__}")
    if isinstance(f, openmm.CustomNonbondedForce):
        print(f"  Energy Expression: {f.getEnergyFunction()}")
        print(f"  Num Particles: {f.getNumParticles()}")
        print(f"  Num Global Params: {f.getNumGlobalParameters()}")
        for k in range(f.getNumGlobalParameters()):
             print(f"    Param {k}: {f.getGlobalParameterName(k)} = {f.getGlobalParameterDefaultValue(k)}")
        
        # Check per-particle params
        if f.getNumParticles() > 0:
             print(f"    Particle 0 params: {f.getParticleParameters(0)}")

    if isinstance(f, openmm.GBSAOBCForce):
        print(f"  GBSAOBCForce found. SoluteDielectric: {f.getSoluteDielectric()}, SolventDielectric: {f.getSolventDielectric()}")


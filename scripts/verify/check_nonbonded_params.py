
import openmm
import openmm.app as app
import sys

# Load prmtop
print("Loading prmtop...")
prmtop = app.AmberPrmtopFile('test/data/6t1h_1151_comp/complex.prmtop')
system = prmtop.createSystem(nonbondedMethod=app.NoCutoff, implicitSolvent=app.OBC2)

print("Inspecting NonbondedForce...")
nb_force = None
cnb_force = None
for f in system.getForces():
    if isinstance(f, openmm.NonbondedForce):
        nb_force = f
    if isinstance(f, openmm.CustomNonbondedForce):
        cnb_force = f

if nb_force:
    print(f"NonbondedForce: {nb_force.getNumParticles()} particles")
    # Check first few particles
    for i in range(min(5, nb_force.getNumParticles())):
        chg, sig, eps = nb_force.getParticleParameters(i)
        print(f"  Particle {i}: Charge={chg}, Sigma={sig}, Epsilon={eps}")

if cnb_force:
    print(f"\nCustomNonbondedForce: {cnb_force.getNumParticles()} particles")
    print(f"Energy function: {cnb_force.getEnergyFunction()}")
    # Check what parameters are
    # Using getParticleParameters(i) returns tuple of params
    # We need to know what they map to.
    print(f"Per-particle parameters: num={cnb_force.getNumPerParticleParameters()}")
    for k in range(cnb_force.getNumPerParticleParameters()):
        print(f"  Param {k}: {cnb_force.getPerParticleParameterName(k)}")
    
    for i in range(min(5, cnb_force.getNumParticles())):
        params = cnb_force.getParticleParameters(i)
        print(f"  Particle {i}: {params}")

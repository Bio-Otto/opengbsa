
import openmm
import sys

print(f"OpenMM Version: {openmm.__version__}")

f = openmm.GBSAOBCForce()
print("Force created: GBSAOBCForce")

# Check for parameter setters (alpha, beta, gamma)
# These are typically not exposed directly as setters on the force object in older versions, 
# but are properties or arguments to constructor in very old ones? 
# Actually in OpenMM 7+, they might be hidden or implicit?
# Let's check dir
print("\nAttributes/Methods:")
for x in dir(f):
    if 'alpha' in x.lower() or 'beta' in x.lower() or 'gamma' in x.lower() or 'kappa' in x.lower():
        print(f"  {x}")

# Getting default values might be hard without particles?
# But let's see if we can get them.
try:
   # usually they are not queryable global params?
   pass
except:
   pass

# Check if we can set them?
# OpenMM documentation says GBSAOBCForce uses OBC2 by default?
# Let's try to find how to switch between OBC1 and OBC2.
# Usually it's via separate parameters.

print("\nKappa check:")
# Salt concentration
if hasattr(f, 'getSolventDielectric'):
    print(f"  getSolventDielectric: {f.getSolventDielectric()}")
if hasattr(f, 'getSoluteDielectric'):
    print(f"  getSoluteDielectric: {f.getSoluteDielectric()}")


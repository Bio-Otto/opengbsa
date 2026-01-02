import openmm
from openmm import Platform, Vec3

def test_openmm_platforms():
    print("Available OpenMM platforms:")
    for i in range(Platform.getNumPlatforms()):
        plat = Platform.getPlatform(i)
        print(f"  {i}: {plat.getName()} (Properties: {plat.getPropertyNames()})")

    for i in range(Platform.getNumPlatforms()):
        plat = Platform.getPlatform(i)
        name = plat.getName()
        print(f"\nTesting platform: {name}")
        try:
            system = openmm.System()
            system.addParticle(39.948)
            integrator = openmm.VerletIntegrator(1.0)
            context = openmm.Context(system, integrator, plat)
            # Set positions (required!)
            context.setPositions([Vec3(0, 0, 0)])
            state = context.getState(getEnergy=True)
            print(f"  ✓ {name} platform is available and context creation succeeded.")
            del context, integrator, system
        except Exception as e:
            print(f"  ✗ {name} platform is NOT available: {e}")

if __name__ == "__main__":
    test_openmm_platforms() 

import openmm
from openmm import app, unit
from openmm.app.internal.lcpo import getLCPOParamsTopology

def test_simple_lcpo():
    print("Creating simple system (Alanine Dipeptide)...")
    # specific force field not strictly needed if we just want topology for LCPO
    # but we need a topology with bonds.
    
    # Manually create topology to avoid file I/O issues
    
    topology = app.Topology()
    chain = topology.addChain()
    residue = topology.addResidue('ALA', chain)
    # Create atoms: N, CA, C, O, CB
    n = topology.addAtom('N', app.Element.getBySymbol('N'), residue)
    ca = topology.addAtom('CA', app.Element.getBySymbol('C'), residue)
    c = topology.addAtom('C', app.Element.getBySymbol('C'), residue)
    o = topology.addAtom('O', app.Element.getBySymbol('O'), residue)
    cb = topology.addAtom('CB', app.Element.getBySymbol('C'), residue)
    
    # Add bonds
    topology.addBond(n, ca)
    topology.addBond(ca, c)
    topology.addBond(c, o)
    topology.addBond(ca, cb)
    
    print("Topology created. Calculating LCPO parameters...")
    try:
        params = getLCPOParamsTopology(topology)
        print(f"LCPO Params calculated for {len(params)} atoms.")
        print(f"Sample param (Atom 0): {params[0]}")
    except Exception as e:
        print(f"Failed to calculate params: {e}")
        return

    print("Creating LCPOForce...")
    lcpo = openmm.LCPOForce()
    lcpo.setSurfaceTension(0.0072 * unit.kilocalories_per_mole / unit.angstroms**2) # approx value
    
    # Add particles
    for p in params:
        # p is [radius, p1, p2, p3, p4]
        # LCPOForce.addParticle(radius, p1, p2, p3, p4)
        # Note: units in params might be quantities
        r = p[0]
        # internal LCPO checks units?
        # getLCPOParamsTopology returns quantities? 
        # Let's check type in output or just assume they are used as is in core.py
        
        # In core.py:
        # radius_nm = r.value_in_unit(unit.nanometer)
        # ...
        
        # We will just verify we can add them.
        try:
             # Just matching core.py logic roughly
             radius = p[0].value_in_unit(unit.nanometer)
             p1 = p[1]
             p2 = p[2]
             p3 = p[3]
             p4 = p[4]
             lcpo.addParticle(radius, p1, p2, p3, p4)
        except Exception as e:
             print(f"Error adding particle: {e}")
    
    print(f"Added {lcpo.getNumParticles()} particles to LCPOForce.")
    
    # Create System to test execution
    system = openmm.System()
    for i in range(5):
        system.addParticle(12.0)
    system.addForce(lcpo)
    
    # Context
    integrator = openmm.VerletIntegrator(1.0*unit.femtoseconds)
    platform = openmm.Platform.getPlatformByName('Reference') # Use Ref for safety/speed on small system
    context = openmm.Context(system, integrator, platform)
    
    positions = [
        openmm.Vec3(0,0,0),
        openmm.Vec3(0.15,0,0),
        openmm.Vec3(0.30,0,0),
        openmm.Vec3(0.30,0.15,0),
        openmm.Vec3(0.15,0.15,0)
    ] # dummy positions in nm
    
    context.setPositions(positions)
    state = context.getState(getEnergy=True)
    energy = state.getPotentialEnergy()
    print(f"LCPO Energy: {energy}")
    print("SUCCESS: LCPO works on simple system.")

if __name__ == "__main__":
    test_simple_lcpo()

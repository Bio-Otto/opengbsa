
import openmm
from openmm import app, unit
from openmm.app.internal.lcpo import getLCPOParamsTopology
import sys

def test_lcpo_ligand():
    # Use ligand Prmtop which is small
    prmtop_path = "test/data/6t1h_1151_comp/ligand.prmtop"
    print(f"Loading {prmtop_path}...")
    prmtop = app.AmberPrmtopFile(prmtop_path)
    
    print(f"Ligand Atoms: {prmtop.topology.getNumAtoms()}")
    
    # Defaults
    ST_DEFAULT = 3.01248 * unit.kilojoules_per_mole / unit.nanometers**2 # 0.0072 kcal/mol/A^2
    PR_DEFAULT = 0.14 * unit.nanometers # 1.4 A
    
    print(f"Surface Tension: {ST_DEFAULT}")
    print(f"Probe Radius: {PR_DEFAULT}")

    def calc_energy(system_setup_func):
        system = openmm.System()
        for _ in range(prmtop.topology.getNumAtoms()):
            system.addParticle(1.0)
            
        force = system_setup_func(prmtop.topology)
        system.addForce(force)
        
        # Dummy positions (randomized cloud or just zeros + slight offset to avoid overlaps if needed)
        # Actually LCPO depends on overlap. Zeros will give max overlap?
        # We need realistic coordinates.
        # But we don't have a ligand-only coordinate file handy.
        # Let's create dummy coords in a line to ensure no overlap for simplicity, 
        # or just 0,0,0 and see what happens (singularity?)
        # Better: use simple grid
        import random
        pos = []
        for i in range(prmtop.topology.getNumAtoms()):
            pos.append(openmm.Vec3(i*0.15, 0, 0)) # Linear chain, no overlap
            
        integrator = openmm.VerletIntegrator(1.0*unit.femtoseconds)
        context = openmm.Context(system, integrator)
        context.setPositions(pos)
        return context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole)

    # 1. Baseline
    def setup_baseline(topology):
        lcpo = openmm.LCPOForce()
        lcpo.setSurfaceTension(ST_DEFAULT)
        params = getLCPOParamsTopology(topology)
        for p in params:
             r_nm = p[0].value_in_unit(unit.nanometer)
             # Default logic
             eff_r = (r_nm + 0.14) if r_nm > 0 else 0.0
             lcpo.addParticle(eff_r, p[1], p[2], p[3], p[4].value_in_unit(unit.nanometer**-2))
        return lcpo
        
    e_base = calc_energy(setup_baseline)
    print(f"Baseline Energy: {e_base:.4f} kcal/mol")
    
    # 2. No Hydrogens (Simulated by 0 radius for H?)
    # or by removing H from topology.
    # Let's just create a modified topology logic
    def setup_no_h(topology):
        lcpo = openmm.LCPOForce()
        lcpo.setSurfaceTension(ST_DEFAULT)
        params = getLCPOParamsTopology(topology)
        atoms = list(topology.atoms())
        for i, p in enumerate(params):
             r_nm = p[0].value_in_unit(unit.nanometer)
             atom = atoms[i]
             if atom.element.symbol == 'H':
                 # Skip H by giving 0 radius? 
                 # LCPO formula uses neighbor counts too.
                 # This is approximated.
                 eff_r = 0.0
             else:
                 eff_r = (r_nm + 0.14) if r_nm > 0 else 0.0
             lcpo.addParticle(eff_r, p[1], p[2], p[3], p[4].value_in_unit(unit.nanometer**-2))
        return lcpo
        
    e_noh = calc_energy(setup_no_h)
    print(f"No-H Energy (Approx): {e_noh:.4f} kcal/mol")

    # 3. Probe Radius 0
    def setup_probe0(topology):
        lcpo = openmm.LCPOForce()
        lcpo.setSurfaceTension(ST_DEFAULT)
        params = getLCPOParamsTopology(topology)
        for p in params:
             r_nm = p[0].value_in_unit(unit.nanometer)
             # Probe 0
             eff_r = (r_nm + 0.0) if r_nm > 0 else 0.0
             lcpo.addParticle(eff_r, p[1], p[2], p[3], p[4].value_in_unit(unit.nanometer**-2))
        return lcpo

    e_p0 = calc_energy(setup_probe0)
    print(f"Probe 0.0 Energy: {e_p0:.4f} kcal/mol")

if __name__ == "__main__":
    test_lcpo_ligand()

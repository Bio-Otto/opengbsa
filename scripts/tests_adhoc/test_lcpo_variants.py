
import openmm
from openmm import app, unit
from openmm.app.internal.lcpo import getLCPOParamsTopology
import sys

def test_lcpo_variants(prmtop_path, gro_path):
    print(f"Loading {prmtop_path} and {gro_path}...")
    prmtop = app.AmberPrmtopFile(prmtop_path)
    if gro_path.endswith('.gro'):
        gro = app.GromACSFile(gro_path)
    else:
        gro = app.PDBFile(gro_path)
    
    # Ensure they match
    if prmtop.topology.getNumAtoms() != gro.topology.getNumAtoms():
        print("Error: Atom counts differ!")
        return

    # 1. Baseline LCPO (With Hydrogens, Default Params)
    print("\n--- Baseline LCPO (With Hydrogens) ---")
    system = openmm.System()
    for _ in range(prmtop.topology.getNumAtoms()):
        system.addParticle(1.0) # Dummy mass
    
    # Core.py settings
    # surface_tension = 3.01248 * unit.kilojoules_per_mole / unit.nanometers**2  (~ 0.0072 kcal/mol/A^2)
    # probe_radius = 1.4 Angstrom
    
    st_kcal = 0.0072 # kcal/mol/A^2
    st_kj = st_kcal * 4.184 * 100 # kJ/mol/nm^2 approx 3.01
    ST_DEFAULT = st_kj * unit.kilojoules_per_mole / unit.nanometers**2
    PR_DEFAULT = 0.14 * unit.nanometers # 1.4 A
    
    print(f"Surface Tension: {ST_DEFAULT}")
    print(f"Probe Radius: {PR_DEFAULT}")

    lcpo = openmm.LCPOForce()
    lcpo.setSurfaceTension(ST_DEFAULT)
    
    params = getLCPOParamsTopology(prmtop.topology)
    for p in params:
        # p[0] is radius (quantity nm)
        # p[4] is p4 (quantity 1/nm^2)
        # others are unitless
        r_nm = p[0].value_in_unit(unit.nanometer)
        eff_r = (r_nm + 0.14) if r_nm > 0 else 0.0
        lcpo.addParticle(eff_r, p[1], p[2], p[3], p[4].value_in_unit(unit.nanometer**-2))
        
    system.addForce(lcpo)
    
    integrator = openmm.VerletIntegrator(1.0*unit.femtoseconds)
    context = openmm.Context(system, integrator)
    context.setPositions(gro.positions)
    
    e_baseline = context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole)
    print(f"Baseline Energy: {e_baseline:.4f} kcal/mol")
    
    # 2. No Hydrogens LCPO
    # Create a new topology without Hydrogens
    print("\n--- No Hydrogens LCPO ---")
    modeller = app.Modeller(prmtop.topology, gro.positions)
    modeller.deleteWater() # Should be dry already but just in case
    
    # Identify Hydrogens
    hs = [atom for atom in modeller.topology.atoms() if atom.element.symbol == 'H']
    modeller.delete(hs)
    
    print(f"Atoms after stripping H: {modeller.topology.getNumAtoms()}")
    
    system_noh = openmm.System()
    for _ in range(modeller.topology.getNumAtoms()):
        system_noh.addParticle(1.0)
        
    lcpo_noh = openmm.LCPOForce()
    lcpo_noh.setSurfaceTension(ST_DEFAULT)
    
    try:
        params_noh = getLCPOParamsTopology(modeller.topology)
        for p in params_noh:
            r_nm = p[0].value_in_unit(unit.nanometer)
            eff_r = (r_nm + 0.14) if r_nm > 0 else 0.0
            lcpo_noh.addParticle(eff_r, p[1], p[2], p[3], p[4].value_in_unit(unit.nanometer**-2))
            
        system_noh.addForce(lcpo_noh)
        context_noh = openmm.Context(system_noh, integrator)
        context_noh.setPositions(modeller.positions)
        e_noh = context_noh.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole)
        print(f"No-Hydrogen Energy: {e_noh:.4f} kcal/mol")
    except Exception as e:
        print(f"Error in No-H calculation: {e}")

    # 3. Parameter Variations (On Full System)
    # Check if Probe Radius matches
    print("\n--- Parameter Sensitivity (Full System) ---")
    
    def test_params(tension_factor=1.0, probe_radius_nm=0.14):
        s = openmm.System()
        for _ in range(prmtop.topology.getNumAtoms()):
            s.addParticle(1.0)
        l = openmm.LCPOForce()
        l.setSurfaceTension(ST_DEFAULT * tension_factor)
        
        ps = getLCPOParamsTopology(prmtop.topology)
        for p in ps:
            r_nm = p[0].value_in_unit(unit.nanometer)
            eff_r = (r_nm + probe_radius_nm) if r_nm > 0 else 0.0
            l.addParticle(eff_r, p[1], p[2], p[3], p[4].value_in_unit(unit.nanometer**-2))
        s.addForce(l)
        c = openmm.Context(s, integrator)
        c.setPositions(gro.positions)
        return c.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole)

    print(f"Energy (Probe=0.0 A): {test_params(probe_radius_nm=0.0):.4f} kcal/mol")
    print(f"Energy (Tension=0.005): {test_params(tension_factor=0.005/0.0072):.4f} kcal/mol")


if __name__ == "__main__":
    prmtop = "test/data/6t1h_1151_comp/complex.prmtop"
    gro = "test/data/6t1h_1151_comp/replica_1/gro/complex_solv_ions.gro" # Contains water, need dry coords
    # Actually, we should use the dry pdb if available, or just load prmtop coords if present.
    # The runner uses xtc + gro. Here let's just use prmtop for structure and some dummy coords if needed
    # Wait, LCPO depends on positions. We need a valid frame.
    # OpenGBSA extracts frames. Let's use 'temp_fixed_for_panda.pdb' from results if it exists
    
    # Hardcoded latest result path for reliability
    pdb_path = "test/results/6t1h_1151_comp/analysis_20260208_230213/temp_fixed_for_panda.pdb"

    if pdb_path and os.path.exists(pdb_path):
        print(f"Using Coords from: {pdb_path}")
        # temp_fixed_for_panda is PDB.
        pdb = app.PDBFile(pdb_path)
        # We need to match it with prmtop. PDB might have lost atom types/bond info needed for LCPO?
        # getLCPOParamsTopology needs bonds and elements. PDB has them.
        # But we want to test on the prmtop topology to match the run.
        # So: Load Prmtop, load PDB positions.
        
        # Check atom count match
        prm = app.AmberPrmtopFile(prmtop)
        if prm.topology.getNumAtoms() == pdb.topology.getNumAtoms():
             test_lcpo_variants(prmtop, pdb_path) # Pass pdb as gro argument for coords
        else:
             print("Atom count mismatch between prmtop and result pdb. Using PDB topology directly for test.")
             # This bypasses the prmtop issue but tests the LCPO code itself
             test_lcpo_variants(pdb_path, pdb_path) # Use pdb for both
            
    else:
        print("Could not find result PDB. Cannot run test.")


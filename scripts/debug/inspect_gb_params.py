
import openmm
from openmm import app, unit
import parmed as pmd
import sys

def inspect_gb_params():
    prmtop_path = 'test/data/6t1h_1151_comp/receptor.prmtop'
    print(f"Loading {prmtop_path}...")
    prmtop = app.AmberPrmtopFile(prmtop_path)
    
    # Create System with OBC2 (igb=5 equivalent)
    system = prmtop.createSystem(implicitSolvent=app.OBC2, 
                                 nonbondedMethod=app.NoCutoff,
                                 constraints=None,
                                 rigidWater=False)
    
    # Find GB Force
    gb_force = None
    for f in system.getForces():
        if isinstance(f, openmm.GBSAOBCForce):
            gb_force = f
            break
            
    if not gb_force:
        print("Error: No GBSAOBCForce found!")
        return

    print(f"GB Force found: {gb_force.__class__.__name__}")
    print(f"Solute Dielectric: {gb_force.getSoluteDielectric()}")
    print(f"Solvent Dielectric: {gb_force.getSolventDielectric()}")
    print(f"Number of Particles: {gb_force.getNumParticles()}")
    
    # Check Exclusions
    # GBSAOBCForce doesn't have getNumExclusions in strict API?
    # It assumes exclusions from NonbondedForce?
    # Actually it has addParticle but usually inherits exceptions/exclusions?
    # Wait, in OpenMM GBSAOBCForce does NOT have exclusions management methods?
    # It uses the system's exclusions?
    # Let's check dir()
    print(f"GB Force methods: {[m for m in dir(gb_force) if 'xclusion' in m]}")
    
    # Check NonbondedForce exclusions
    nb_force = [f for f in system.getForces() if isinstance(f, openmm.NonbondedForce)][0]
    print(f"Nonbonded Force Exclusions: {nb_force.getNumExceptions()}")
    
    # Load with ParmEd to check intrinsic radii
    
    # Load with ParmEd to check intrinsic radii
    struct = pmd.load_file(prmtop_path)
    
    print("\nChecking first 10 atoms...")
    print(f"{'Atom':<10} {'Elem':<5} {'Prmtop Rad':<12} {'Prmtop Screen':<12} {'OpenMM Rad':<12} {'OpenMM Scale':<12} {'Diff Rad'}")
    print("-" * 80)
    
    # Check first few atoms
    for i in range(10):
        # OpenMM Params
        q, rad_nm, scale = gb_force.getParticleParameters(i)
        rad_ang = rad_nm.value_in_unit(unit.angstroms)
        
        # ParmEd Params
        atom = struct.atoms[i]
        p_rad = atom.solvent_radius # Angstroms
        p_screen = atom.screen
        
        diff = abs(rad_ang - p_rad)
        
        print(f"{atom.name:<10} {atom.element_name:<5} {p_rad:<12.4f} {p_screen:<12.4f} {rad_ang:<12.4f} {scale:<12.4f} {diff:.4f}")

    # Check summary
    total_diff = 0
    max_diff = 0
    for i in range(system.getNumParticles()):
         q, rad_nm, scale = gb_force.getParticleParameters(i)
         rad_ang = rad_nm.value_in_unit(unit.angstroms)
         atom = struct.atoms[i]
         diff = abs(rad_ang - atom.solvent_radius)
         total_diff += diff
         if diff > max_diff:
             max_diff = diff
             
    print(f"\nSummary:")
    print(f"Total Radii Diff (Sum): {total_diff:.4f}")
    print(f"Max Radii Diff: {max_diff:.4f}")
    
    if max_diff < 0.001:
        print("✓ OpenMM uses Prmtop radii correctly.")
    else:
        print("❌ OpenMM Radii Mismatch! System creation is ignoring Prmtop radii?")

if __name__ == "__main__":
    inspect_gb_params()

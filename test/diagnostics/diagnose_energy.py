
import sys
import os
import yaml
import numpy as np
import mdtraj as md
from mmgbsa.runner import MMGBSARunner
from mmgbsa.core import FixedEnhancedTrueForceFieldMMGBSA
from openmm import unit, app

def diagnose_frame(config_file):
    print(f"Diagnosing configuration: {config_file}")
    
    # Load config
    with open(config_file) as f:
        config = yaml.safe_load(f)
        
    runner = MMGBSARunner(config_file)
    input_files = runner.config['input_files']
    
    complex_pdb = input_files['complex_pdb']
    xtc_file = input_files['trajectory']
    ligand_mol = input_files['ligand_mol']
    ligand_pdb = input_files['ligand_pdb']
    
    print(f"Loading trajectory: {xtc_file}")
    traj = md.load(xtc_file, top=complex_pdb)
    print(f"Loaded {len(traj)} frames.")
    
    # We only care about frame 0 for diagnosis
    frame = traj[0]
    


    # Setup calculator directly
    calculator = FixedEnhancedTrueForceFieldMMGBSA(
        temperature=float(config['analysis_settings'].get('temperature', 300)) * unit.kelvin,
        gb_model=config['analysis_settings'].get('gb_model', 'OBC2'),
        salt_concentration=float(config['analysis_settings'].get('salt_concentration', 0.15)) * unit.molar
    )
    
    # Setup system (this uses the same logic as the main run)
    print("Setting up system...")
    platform, props = calculator.setup_optimized_platform()
    
    # Load complex topology
    pdb = app.PDBFile(complex_pdb)
    complex_top = pdb.topology
    ligand_resname = calculator.find_ligand_resname(complex_top)
    print(f"Ligand resname: {ligand_resname}")
    
    ligand_indices = calculator.get_ligand_indices(complex_top, ligand_resname)
    protein_indices = calculator.get_protein_indices(complex_top, ligand_resname)
    
    # Parameterize
    ligand_system, ligand_top, ligand_mol_obj = calculator.parameterize_ligand_openff(ligand_mol)
    protein_system, protein_top, protein_pos = calculator.parameterize_protein_amber(complex_pdb, ligand_resname)
    complex_system = calculator.create_combined_system(protein_system, ligand_system)
    
    # Context
    complex_integrator = openmm.LangevinMiddleIntegrator(300*unit.kelvin, 1/unit.picosecond, 0.001*unit.picosecond)
    complex_context = openmm.Context(complex_system, complex_integrator, platform, props)
    
    # Set positions from trajectory Frame 0
    xyz = frame.xyz[0] * unit.nanometer
    complex_context.setPositions(xyz)
    
    state = complex_context.getState(getEnergy=True, getForces=True)
    pot_e = state.getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole)
    print(f"\nTotal Potential Energy (Frame 0): {pot_e:.2f} kcal/mol")
    
    # Decompose components
    # Nonbonded is usually Group 0 or 1 depending on implementation
    # In core.py: 
    # Group 1 = GBSAOBCForce (GB)
    # Group 0 = NonbondedForce (vdW + Elec)
    
    e_nb = complex_context.getState(getEnergy=True, groups=1).getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole) # GB? No check core.py
    # core.py lines 1607+:
    # e_nb = groups=1 ?? Wait, standard NonbondedForce is usually group 0.
    # Let's check groups in system.
    
    forces = {f.getName(): f for f in complex_system.getForces()}
    print("\nForce Groups:")
    for f in complex_system.getForces():
        print(f"  {f.getName()}: Group {f.getForceGroup()}")
        
    # Recalculate based on specific groups
    for group_idx in range(32):
        try:
            e = complex_context.getState(getEnergy=True, groups=(1<<group_idx)).getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole)
            if abs(e) > 0.001:
                print(f"Group {group_idx} Energy: {e:.2f} kcal/mol")
        except:
            pass
            
    # Check Geometry
    lig_xyz = xyz[ligand_indices]
    prot_xyz = xyz[protein_indices]
    
    lig_com = np.mean(lig_xyz.value_in_unit(unit.nanometer), axis=0)
    prot_com = np.mean(prot_xyz.value_in_unit(unit.nanometer), axis=0)
    
    dist = np.linalg.norm(lig_com - prot_com)
    print(f"\nDistance between Protein COM and Ligand COM: {dist:.2f} nm")
    
    # Check for Clashes
    # Simple check: min distance between any protein atom and any ligand atom
    print("Checking min atomic distance...")
    min_dist = 9999.0
    
    # This is slow O(N*M), let's do a quick naive implementation or use mdtraj
    t_frame = traj[0]
    pairs = t_frame.top.select_pairs(f"resid {ligand_resname}", "protein")
    if len(pairs) > 0:
        dists = md.compute_distances(t_frame, pairs)
        min_dist = np.min(dists)
        print(f"Minimum Protein-Ligand Atom Distance: {min_dist:.4f} nm")
        
        if min_dist < 0.1: # less than 1 Angstrom
            print("⚠️  CRITICAL: Severe atomic clash detected! Atoms are overlapping.")
        elif min_dist < 0.2:
            print("⚠️  Warning: Very close contact detected.")
        else:
            print("✅ Geometry looks reasonable (no obvious clashes).")
    else:
        print("Could not compute distances (selection failed?)")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python diagnos_energy.py config.yaml")
        sys.exit(1)
        
    diagnose_frame(sys.argv[1])

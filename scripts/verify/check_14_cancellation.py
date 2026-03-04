
import parmed as pmd
from openmm import app, openmm, unit
import mdtraj as md
import numpy as np

def check_14_cancellation():
    print("Checking 1-4 Energy Cancellation...")
    prmtop = "test/data/6t1h_6466_comp/complex.prmtop"
    xtc = "test/data/6t1h_6466_comp/md_complex_prod.xtc"
    solv_top = "test/data/6t1h_6466_comp/replica_1/gro/complex_solv_ions.gro"
    
    struct = pmd.load_file(prmtop)
    
    # Split using ParmEd (Logic from Core.py)
    # Note: Core.py uses structure.copy() and strip()
    ligand_resname = "LIG"
    
    print("Splitting structures...")
    complex_struct = struct
    receptor_struct = struct.copy(cls=pmd.Structure)
    receptor_struct.strip(f":{ligand_resname}")
    ligand_struct = struct.copy(cls=pmd.Structure)
    ligand_struct.strip(f"!:{ligand_resname}")
    
    # Load Weights (Positions)
    print("Loading positions...")
    traj = md.load_frame(xtc, 0, top=solv_top)
    all_pos = traj.xyz[0] * unit.nanometer
    
    # Indices for subsets
    # Match by index? strip() removes atoms, so indices change.
    # We need to map original indices to new indices.
    # Or just use the positions directly from the subsets?
    # No, we need consistent positions.
    
    # Get mask for receptor/ligand in original
    lig_mask = [a.idx for a in struct.atoms if a.residue.name == ligand_resname]
    rec_mask = [a.idx for a in struct.atoms if a.residue.name != ligand_resname]
    
    complex_pos = all_pos[:len(struct.atoms)] # First N atoms
    rec_pos = complex_pos[rec_mask]
    lig_pos = complex_pos[lig_mask]
    
    # Helper to calc 1-4 E
    def calc_e14(structure, positions, label):
        sys = structure.createSystem(nonbondedMethod=app.NoCutoff)
        nb = [f for f in sys.getForces() if isinstance(f, openmm.NonbondedForce)][0]
        
        # Zero non-1-4s
        for i in range(nb.getNumParticles()):
            nb.setParticleParameters(i, 0.0, 0.1, 0.0)
            
        ctx = openmm.Context(sys, openmm.VerletIntegrator(1.0*unit.femtoseconds))
        ctx.setPositions(positions)
        e = ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole)
        print(f"  {label} 1-4 E: {e:.4f} kcal/mol")
        return e

    e_complex = calc_e14(complex_struct, complex_pos, "Complex")
    e_receptor = calc_e14(receptor_struct, rec_pos, "Receptor")
    e_ligand = calc_e14(ligand_struct, lig_pos, "Ligand")
    
    delta = e_complex - e_receptor - e_ligand
    print(f"\nCancellation Delta: {delta:.4f} kcal/mol")
    
    if abs(delta) > 0.1:
        print("FAIL: Significant cancellation error!")
        print("This explains why VDW/ELE energies are wrong.")
    else:
        print("PASS: 1-4 energies cancel perfectly.")

if __name__ == "__main__":
    check_14_cancellation()

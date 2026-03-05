
import openmm
import openmm.app as app
import openmm.unit as unit
import numpy as np
import pandas as pd
import mdtraj as md
import logging

log = logging.getLogger(__name__)

class EnergyDecomposer:
    """
    Efficiently calculates per-residue interaction energies (VdW + Electrostatics)
    between a Ligand and Protein using Numpy vectorization.
    """
    
    def __init__(self, system, topology, protein_indices, ligand_indices, parmed_structure=None):
        self.system = system
        self.topology = topology
        self.protein_indices = protein_indices
        self.ligand_indices = ligand_indices
        self.parmed_structure = parmed_structure
        
        # Extract Parameters (Charge, Sigma, Epsilon)
        self.atom_params = self._extract_parameters()
        
        # Pre-filter protein residues (map residue index to atom indices)
        self.residue_map = self._map_residues()
        
    def _extract_parameters(self):
        """Extract q, sigma, epsilon for all atoms from ParmEd or NonbondedForce"""
        params = {}
        
        # Priority: ParmEd Structure (Robust for Charmm/Gromacs conversions)
        if self.parmed_structure:
            try:
                log.info("FastDecomp: Using ParmEd structure for parameters")
                for i, atom in enumerate(self.parmed_structure.atoms):
                    # Convert Rmin/2 to Sigma (nm)
                    # Sigma = (Rmin/2) / (0.5 * 2^(1/6)) = (Rmin/2) * 1.7818
                    sigma = atom.rmin * 1.781797697 * 0.1 # Angstrom to nm
                    epsilon = atom.epsilon # kcal/mol
                    charge = atom.charge
                    
                    params[i] = {
                        'q': charge,
                        'sigma': sigma,
                        'epsilon': epsilon
                    }
                return params
            except Exception as e:
                log.warning(f"FastDecomp: Failed to extract from ParmEd: {e}")
        
        # Fallback: NonbondedForce
        # Find NonbondedForce
        nb_force = None
        for f in self.system.getForces():
            if isinstance(f, openmm.NonbondedForce):
                nb_force = f
                break
        
        if not nb_force:
            log.warning("No NonbondedForce found for decomposition!")
            return None
            
        for i in range(nb_force.getNumParticles()):
            charge, sigma, epsilon = nb_force.getParticleParameters(i)
            params[i] = {
                'q': charge.value_in_unit(unit.elementary_charge),
                'sigma': sigma.value_in_unit(unit.nanometer),
                'epsilon': epsilon.value_in_unit(unit.kilocalories_per_mole)
            }
        return params

    def _map_residues(self):
        """Map residue objects to their atom indices"""
        res_map = {}
        for res in self.topology.residues():
            # Check if residue is in protein indices
            atom_indices = [a.index for a in res.atoms()]
            if any(idx in self.protein_indices for idx in atom_indices):
                # Filter only protein atoms
                clean_indices = [idx for idx in atom_indices if idx in self.protein_indices]
                if clean_indices:
                    res_map[res] = clean_indices
        return res_map

    def calculate_decomposition(self, positions):
        """
        Calculate interaction energy for current frame positions.
        positions: (N, 3) numpy array in nanometers
        """
        if self.atom_params is None:
            return {}
            
        results = {}
        
        # Constants
        # COULOMB_CONST = 332.0522 kcal*A / ...
        # But r is in nm. So we use 33.20522
        COULOMB_CONST = 33.20522
        
        pos_ligand = positions[self.ligand_indices]
        
        # Pre-calc Ligand Params
        lig_q = np.array([self.atom_params[i]['q'] for i in self.ligand_indices])
        lig_sig = np.array([self.atom_params[i]['sigma'] for i in self.ligand_indices])
        lig_eps = np.array([self.atom_params[i]['epsilon'] for i in self.ligand_indices])
        
        # Vectorized Loop over Residues
        for res, atom_idxs in self.residue_map.items():
            pos_res = positions[atom_idxs]
            
            # Params for Residue
            res_q = np.array([self.atom_params[i]['q'] for i in atom_idxs])
            res_sig = np.array([self.atom_params[i]['sigma'] for i in atom_idxs])
            res_eps = np.array([self.atom_params[i]['epsilon'] for i in atom_idxs])
            
            # Pairwise Distances (N_res x N_lig)
            # Broadcasting: (N_res, 1, 3) - (1, N_lig, 3)
            delta = pos_res[:, np.newaxis, :] - pos_ligand[np.newaxis, :, :]
            r2 = np.sum(delta**2, axis=2)
            r = np.sqrt(r2)
            
            # Avoid division by zero
            r[r < 0.05] = 0.05 
            
            # Electrostatics: k * q1 * q2 / r
            qq = np.outer(res_q, lig_q)
            e_ele = COULOMB_CONST * qq / r
            
            # LJ: 4*eps * ((sig/r)^12 - (sig/r)^6)
            sig_ij = 0.5 * (res_sig[:, np.newaxis] + lig_sig[np.newaxis, :])
            eps_ij = np.sqrt(np.outer(res_eps, lig_eps))
            
            sr6 = (sig_ij / r)**6
            sr12 = sr6**2
            e_vdw = 4.0 * eps_ij * (sr12 - sr6)
            
            total_ele = np.sum(e_ele)
            total_vdw = np.sum(e_vdw)
            
            # Store if significant
            if abs(total_ele + total_vdw) > 0.01:
                # Use PDB-style naming: RES_NUM_CHAIN (e.g. GLU_11_A)
                # Ensure we handle insertion codes if present in id, though OpenMM usually stores it in id
                res_key = f"{res.name}_{res.id}_{res.chain.id}"
                results[res_key] = {
                    'electrostatic': total_ele,
                    'vdw': total_vdw,
                    'total': total_ele + total_vdw
                }
                
        return results

    def analyze_trajectory(self, trajectory):
        """Run decomposition on all frames"""
        logging.info("Starting Per-Residue Decomposition...")
        data = []
        for i in range(len(trajectory)):
            pos = trajectory.xyz[i] # nm
            frame_res = self.calculate_decomposition(pos)
            data.append(frame_res)
            
        return data

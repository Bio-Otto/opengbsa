"""
Topology and structure management for MM/GBSA calculations.
Handles system building, atom selection, and structure manipulation.
"""

import os
import numpy as np
import mdtraj as md
import parmed as pmd
import openmm as mm
import openmm.app as app


class TopologyManager:
    """
    Manages system topology, structure building, and atom selection.
    """
    
    def __init__(self, verbose=False):
        """
        Initialize topology manager.
        
        Parameters:
        -----------
        verbose : bool
            Enable verbose logging
        """
        self.verbose = verbose
    
    def find_ligand_resname(self, complex_pdb, protein_residues=None):
        """
        Identify ligand residue name by process of elimination.
        
        Parameters:
        -----------
        complex_pdb : str
            Path to complex PDB file
        protein_residues : list, optional
            Known protein residue names
            
        Returns:
        --------
        str : Likely ligand residue name or None
        """
        if not os.path.exists(complex_pdb):
            return None
        
        try:
            structure = pmd.load_file(complex_pdb)
            all_residues = set(res.name for res in structure.residues)
            
            # Common protein residues
            common_protein = {
                'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY',
                'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER',
                'THR', 'TRP', 'TYR', 'VAL', 'HIP', 'HID', 'HIE'
            }
            
            candidates = all_residues - common_protein
            
            if candidates:
                ligand_name = sorted(list(candidates))[0]
                if self.verbose:
                    print(f"  Detected ligand residue: {ligand_name}")
                return ligand_name
            
            return None
        except Exception as e:
            if self.verbose:
                print(f"⚠ Error identifying ligand: {e}")
            return None
    
    def get_selection_indices(self, topology, selection_string):
        """
        Get atom indices matching MDTraj selection string.
        
        Parameters:
        -----------
        topology : mdtraj.Topology
            MDTraj topology
        selection_string : str
            MDTraj selection string (e.g., 'protein', 'chainid 0')
            
        Returns:
        --------
        np.ndarray : Indices of selected atoms or None
        """
        if not topology or not selection_string:
            return None
        
        try:
            indices = topology.select(selection_string)
            if self.verbose and indices is not None:
                print(f"  Selection '{selection_string}': {len(indices)} atoms")
            return indices
        except Exception as e:
            if self.verbose:
                print(f"⚠ Selection failed: '{selection_string}' - {e}")
            return None
    
    def get_protein_indices(self, topology):
        """
        Get indices of protein atoms.
        
        Parameters:
        -----------
        topology : mdtraj.Topology
            MDTraj topology
            
        Returns:
        --------
        np.ndarray : Protein atom indices
        """
        return self.get_selection_indices(topology, 'protein')
    
    def get_ligand_indices(self, topology, ligand_resname=None, ligand_selection=None):
        """
        Get indices of ligand atoms.
        
        Parameters:
        -----------
        topology : mdtraj.Topology
            MDTraj topology
        ligand_resname : str, optional
            Ligand residue name
        ligand_selection : str, optional
            MDTraj selection string for ligand
            
        Returns:
        --------
        np.ndarray : Ligand atom indices
        """
        if ligand_selection:
            return self.get_selection_indices(topology, ligand_selection)
        elif ligand_resname:
            return self.get_selection_indices(topology, f'resname {ligand_resname}')
        else:
            return None
    
    def build_complex_system(self, complex_pdb, forcefield=None, implicit_solvent=None):
        """
        Build OpenMM system from complex PDB.
        
        Parameters:
        -----------
        complex_pdb : str
            Path to complex PDB
        forcefield : str, optional
            Forcefield name (e.g., 'amber99sb', 'charmm36')
        implicit_solvent : str, optional
            Implicit solvent model (e.g., 'OBC2', 'GBn')
            
        Returns:
        --------
        tuple : (system, topology, positions) or (None, None, None) on error
        """
        if not os.path.exists(complex_pdb):
            if self.verbose:
                print(f"✗ Complex PDB not found: {complex_pdb}")
            return None, None, None
        
        try:
            # Load PDB
            pdb = app.PDBFile(complex_pdb)
            topology = pdb.topology
            positions = pdb.positions
            
            # Select forcefield
            if not forcefield:
                forcefield = 'amber99sb-ildn'
            
            ff = app.ForceField(f'{forcefield}.xml')
            
            # Create system
            system = ff.createSystem(
                topology,
                nonbondedMethod=app.NoCutoff,
                implicitSolvent=None,
                constraints=app.HBonds
            )
            
            if self.verbose:
                print(f"✓ Built system: {topology.getNumAtoms()} atoms")
            
            return system, topology, positions
            
        except Exception as e:
            if self.verbose:
                print(f"✗ Error building system: {e}")
            return None, None, None
    
    def load_prmtop_complex(self, prmtop_path, gro_or_pdb_path=None, xtc_path=None, 
                           solvated_topology=None, stride=1):
        """
        Load Amber PRMTOP complex topology with trajectory.
        
        Parameters:
        -----------
        prmtop_path : str
            Path to PRMTOP file
        gro_or_pdb_path : str, optional
            Path to initial structure (GRO or PDB)
        xtc_path : str, optional
            Path to XTC trajectory
        solvated_topology : str, optional
            Path to solvated topology for loading XTC
        stride : int
            Stride for loading trajectory
            
        Returns:
        --------
        dict : {'topology': Topology, 'structure': ParmEd Structure, 'traj': MDTraj Trajectory}
        """
        if not os.path.exists(prmtop_path):
            if self.verbose:
                print(f"✗ PRMTOP not found: {prmtop_path}")
            return None
        
        try:
            # Load PRMTOP
            structure = pmd.load_file(prmtop_path)
            
            result = {
                'structure': structure,
                'parmed_structure': structure
            }
            
            # Load trajectory if provided
            if xtc_path and os.path.exists(xtc_path):
                top_file = solvated_topology if (solvated_topology and os.path.exists(solvated_topology)) else prmtop_path
                try:
                    traj = md.load(xtc_path, top=top_file, stride=stride)
                    result['trajectory'] = traj
                except Exception as e:
                    if self.verbose:
                        print(f"⚠ Could not load trajectory: {e}")
            
            if self.verbose:
                print(f"✓ Loaded PRMTOP: {len(structure.atoms)} atoms")
            
            return result
            
        except Exception as e:
            if self.verbose:
                print(f"✗ Error loading PRMTOP: {e}")
            return None

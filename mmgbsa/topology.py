import openmm.app as app
import openmm.unit as unit
from pathlib import Path
from .inputs import EngineMode

class TopologyLoader:
    """
    Factory for creating OpenMM Systems from various topology file formats.
    """
    
    @staticmethod
    def load_system(topology_file: str, mode: EngineMode, **kwargs):
        """
        Loads a topology file and creates an OpenMM System.
        
        Args:
            topology_file: Path to the topology file.
            mode: The EngineMode detected for this file.
            **kwargs: Additional arguments for createSystem (e.g., nonbondedMethod, constraints).
        
        Returns:
            tuple: (openmm.System, openmm.Topology, openmm.unit.Quantity or None)
        """
        path = str(topology_file)
        
        # Default system options (can be overridden by kwargs)
        sys_options = {
            'nonbondedMethod': app.NoCutoff,
            'constraints': None,
            'implicitSolvent': None # strictly vacuum base, add GBSA later
        }
        sys_options.update(kwargs)
        
        if mode == EngineMode.AMBER:
            return TopologyLoader._load_amber(path, sys_options)
        elif mode == EngineMode.GROMACS:
            return TopologyLoader._load_gromacs(path, sys_options)
        elif mode == EngineMode.CHARMM:
            return TopologyLoader._load_charmm(path, sys_options)
        elif mode == EngineMode.OPENMM:
            return TopologyLoader._load_openmm_xml(path, sys_options)
        elif mode == EngineMode.GENERIC:
            return TopologyLoader._load_generic(path, sys_options)
        else:
            raise NotImplementedError(f"System loading for mode {mode} not fully implemented yet.")

    @staticmethod
    def _load_openmm_xml(path, options):
        """Loads OpenMM Serialized System XML."""
        with open(path, 'r') as f:
            xml_content = f.read()
        
        # Determine if it's a System or something else?
        # Typically one saves System via XmlSerializer.
        import openmm as mm
        system = mm.XmlSerializer.deserialize(xml_content)
        
        # We also need a topology for atom selections.
        # XML system doesn't store Topology/Residues names usually.
        # User MUST provide a PDB via 'complex_pdb' logic in core, 
        # but TopologyLoader needs to return a topology.
        # If .xml is provided as 'topology', we can't extract topology easily.
        # Workaround: Warn user that they must rely on PDB for topology logic or provide pdb as companion.
        # But looking at core.py, complex_pdb is often passed to load_system.
        # If complex_pdb is .xml, we have no topology.
        # We return None for topology, but core might crash.
        # Let's check provided companion files if possible? No.
        # Suggestion: Require .pdb for topology if using .xml for system.
        
        # For now, return None for topology and let caller handle it (core usage of TopologyLoader)
        # Actually core line 2097: _, complex_top, _ = TopologyLoader.load_system(complex_pdb, mode)
        # So we NEED a topology.
        # If the user uses .xml, they are likely stuck unless we find a PDB.
        
        # Let's assume the user matches xml name with pdb?
        # Or maybe raise error that XML topology is not sufficient for identification?
        
        # Let's look for a pdb with valid names in same dir?
        p = Path(path)
        candidates = list(p.parent.glob(f"{p.stem}.pdb"))
        topology = None
        if candidates:
             pdb = app.PDBFile(str(candidates[0]))
             topology = pdb.topology
        else:
             # Create dummy topology? No.
             print("WARNING: Loading XML System but could not find matching .pdb for Topology information.")
             
        return system, topology, None

    @staticmethod
    def _load_amber(path, options):
        """Loads Amber prmtop."""
        prmtop = app.AmberPrmtopFile(path)
        system = prmtop.createSystem(**options)
        return system, prmtop.topology, None

    @staticmethod
    def _load_gromacs(path, options):
        """Loads Gromacs top/tpr."""
        # OpenMM GromacsTopFile handles .top, GromacsGroFile handles .gro
        # .tpr support usually requires conversion or specific tools.
        # OpenMM's GromacsTopFile is strictly for .top
        # If user provides .tpr, we might need 'gmx_MMPBSA' style conversion or Parmed.
        # For now, let's assume .top for direct OpenMM support or raise warning.
        if path.endswith('.tpr'):
             # Placeholder for TPR support (ParmEd?)
             import parmed as pmd
             struct = pmd.load_file(path)
             system = struct.createSystem(**options)
             return system, struct.topology, struct.positions
        else:
             top = app.GromacsTopFile(path)
             # Filter options for GromacsTopFile.createSystem
             gmx_options = {k: v for k, v in options.items() if k in ['nonbondedMethod', 'nonbondedCutoff', 'constraints', 'rigidWater', 'removeCMMotion', 'hydrogenMass']}
             system = top.createSystem(**gmx_options)
             return system, top.topology, None

    @staticmethod
    def _load_charmm(path, options):
        """Loads Charmm PSF."""
        psf = app.CharmmPsfFile(path)
        # Charmm needs parameter files loaded into a ForceField logic usually,
        # or PsfFile can create system if params are provided? 
        # app.CharmmPsfFile.createSystem needs params.
        # This is more complex. Sticking to basic stub.
        raise NotImplementedError("Charmm loading requires parameter files.")

    @staticmethod
    def _load_generic(path, options):
        """
        Loads a generic PDB file and constructs a system using standard ForceFields.
        Default: Amber14 + OpenFF/GAFF for ligands if needed (handled by core usually, 
        but here we build the protein part).
        """
        pdb = app.PDBFile(path)
        
        # Load standard ForceField
        # TODO: Make this configurable via options?
        ff_name = options.get('protein_forcefield', 'amber14-all.xml')
        solvent_ff = options.get('solvent_forcefield', 'amber14/tip3p.xml')
        
        try:
            ff = app.ForceField(ff_name, solvent_ff)
            system = ff.createSystem(pdb.topology, **options)
            return system, pdb.topology, pdb.positions
        except Exception as e:
            # Often fails due to missing residues (ligands)
            # In Generic Mode, if PDB has ligands, standard ForceField fails.
            # We might need to return just the PDB topology and let core handle parameterization?
            # But TopologyLoader promises a System.
            # If implementation fails, we raise understandable error.
            raise ValueError(f"Failed to build system from PDB using {ff_name}: {e}. "
                             "Generic PDB mode requires standard residues or explicitly provided template generators.")

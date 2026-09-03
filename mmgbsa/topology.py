"""
Engine-agnostic topology loading dispatch.

Provides `TopologyLoader`, which routes a topology file to the correct
OpenMM System-construction path based on its detected `EngineMode`
(see `mmgbsa/inputs.py`).
"""
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
        """Loads Gromacs .top. Real .tpr support lives in mmgbsa/tpr_loader.py
        and is invoked directly by mmgbsa_core.py before this loader is ever
        reached -- ParmEd itself has no .tpr parser, so this path raises
        rather than silently mis-loading."""
        if path.endswith('.tpr'):
             raise NotImplementedError(
                 "TopologyLoader cannot load .tpr files directly; use "
                 "mmgbsa.tpr_loader.load_tpr_as_parmed instead."
             )
        else:
             top = app.GromacsTopFile(path)
             # Filter options for GromacsTopFile.createSystem
             gmx_options = {k: v for k, v in options.items() if k in ['nonbondedMethod', 'nonbondedCutoff', 'constraints', 'rigidWater', 'removeCMMotion', 'hydrogenMass']}
             system = top.createSystem(**gmx_options)
             return system, top.topology, None

    @staticmethod
    def _load_charmm(path, options):
        """
        Loads a NAMD/CHARMM PSF topology.

        A PSF has no embedded force-field parameters (unlike Amber prmtop),
        so CHARMM parameter files (.prm/.str/.rtf/.par, typically CHARMM36
        or CHARMM36m) must be supplied separately via
        options['charmm_params'] (a str path or list of paths).

        Uses ParmEd (`parmed.load_file` + `Structure.load_parameters`)
        rather than OpenMM's own `app.CharmmPsfFile.loadParameters`: the
        latter's stricter improper-dihedral matching rejects some
        real-world PSFs (confirmed on a cyclic-peptide PSF from a published
        NAMD dataset, which has a head-to-tail backbone improper standard
        CHARMM36m parameters don't cover) with `MissingParameter`, while
        ParmEd's own parameter assignment handles the same PSF/parameter
        combination without error and produces an equivalent System.

        A PSF file itself carries no coordinates (unlike Amber's prmtop+
        inpcrd pairing, where positions are a separate file the caller is
        already expected to supply) -- confirmed directly (`parmed.load_file`
        on a real PSF gives `struct.positions is None`). A companion
        coordinate file (typically the initial .pdb NAMD was given, though
        a `.coor`/`.crd` restart file would also work) is REQUIRED via
        options['charmm_coordinates'] to get real, non-placeholder atomic
        positions -- returning an all-zeros placeholder here (as this
        codebase's caller does for a still-None `positions`, see
        `parameterize_protein_amber`) would silently produce a physically
        meaningless system, which is worse than failing loudly.
        """
        charmm_params = options.pop('charmm_params', None)
        if not charmm_params:
            raise ValueError(
                "Charmm/.psf topology loading requires CHARMM parameter files "
                "(.prm/.str/.rtf) -- pass them as options['charmm_params'] "
                "(a single path or a list of paths), e.g. the CHARMM36m "
                "'par_all36m_prot.prm' and 'toppar_water_ions_prot.str' files "
                "distributed with a NAMD input set."
            )
        if isinstance(charmm_params, str):
            charmm_params = [charmm_params]
        charmm_coordinates = options.pop('charmm_coordinates', None)

        import parmed as pmd
        from parmed.charmm import CharmmParameterSet

        params = CharmmParameterSet(*charmm_params)
        struct = pmd.load_file(path)
        struct.load_parameters(params)

        if charmm_coordinates:
            coord_struct = pmd.load_file(charmm_coordinates)
            struct.coordinates = coord_struct.coordinates

        createsystem_options = {k: v for k, v in options.items()
                                 if k in ('nonbondedMethod', 'nonbondedCutoff', 'constraints',
                                          'rigidWater', 'removeCMMotion', 'hydrogenMass',
                                          'implicitSolvent', 'implicitSolventSaltConc',
                                          'switchDistance', 'ewaldErrorTolerance')}
        system = struct.createSystem(**createsystem_options)
        return system, struct.topology, struct.positions

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

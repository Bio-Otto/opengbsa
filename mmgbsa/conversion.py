"""
GROMACS-to-Amber format conversion via ParmEd.

Provides `GromacsPreprocessor`, which loads a GROMACS `.top`/`.gro` system
and writes out equivalent Amber `.prmtop`/`.inpcrd` files (applying the
configured GB radii set along the way), so downstream code can treat a
GROMACS-origin system uniformly with native Amber input.
"""
import parmed as pmd
import logging
from pathlib import Path
import os
import sys

log = logging.getLogger(__name__)

class GromacsPreprocessor:
    """
    Handles automatic conversion of Gromacs files (.top, .gro) to Amber format (.prmtop)
    to enable Native Amber Mode processing in MM/GBSA.
    """
    
    # Standard Amber GB-model-to-radii-set pairing (see e.g. Amber Reference
    # Manual / MMPBSA.py docs): HCT/GBn pair with 'mbondi', OBC1/OBC2 with
    # 'mbondi2', GBn2 with 'mbondi3'. Kept in sync with the identical mapping
    # in mmgbsa_core.py's StructureManager._GB_MODEL_RADII_SET.
    _GB_MODEL_RADII_SET = {
        'HCT': 'mbondi',
        'OBC1': 'mbondi2',
        'OBC2': 'mbondi2',
        'GBn': 'mbondi',
        'GBn2': 'mbondi3',
    }

    @staticmethod
    def convert_to_amber(topology_file, coordinate_file, output_dir=None, gb_model='OBC2'):
        """
        Converts Gromacs topology to Dry Amber/Parm7 format.

        GROMACS topologies carry no intrinsic GB radii/screen values at all
        (that's an Amber-ecosystem convention, not part of the GROMACS force
        field) -- ParmEd's `GromacsTopologyFile` loader leaves every atom's
        `solvent_radius`/`screen` at 0.0. Saving directly to .prmtop without
        assigning a real radii set means every downstream consumer of these
        converted files silently falls back to a crude, atom-type-blind
        element-only radius table, producing GB energies that are NOT
        comparable to genuine Amber MMPBSA.py results (confirmed: this
        caused a systematic ~20 kcal/mol discrepancy against Amber reference
        results on real OXA-MD benchmark systems before this fix). The GB
        radii set matching `gb_model` is applied via ParmEd's own
        (Amber-validated) `changeRadii` action before any file is saved, so
        every one of the four output .prmtop files carries correct radii.

        Parameters
        ----------
        gb_model : str
            Configured GB model (OBC1/OBC2/HCT/GBn/GBn2), used to select the
            matching radii set (see `_GB_MODEL_RADII_SET`).

        Returns:
            dict: Paths to generated files:
                - 'topology': complex_dry.prmtop (for simulation)
                - 'solvated': complex.prmtop (for trajectory mapping)
                - 'receptor': receptor.prmtop
                - 'ligand': ligand.prmtop
        """
        topology_file = Path(topology_file).resolve()
        coordinate_file = Path(coordinate_file).resolve()
        
        # The converted files go to output_dir (the run folder), never next to the
        # input; without one, to a new temporary folder
        if output_dir is None:
            import tempfile
            output_dir = Path(tempfile.mkdtemp(prefix='opengbsa_gromacs_'))
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
            
        stem = topology_file.stem
        
        # Output paths
        out_solvated = output_dir / f"{stem}_solvated.prmtop"
        out_dry = output_dir / f"{stem}_dry.prmtop"
        out_receptor = output_dir / "receptor.prmtop"
        out_ligand = output_dir / "ligand.prmtop"
        
        log.info(f"Auto-converting Gromacs files: {topology_file.name} -> Amber format")
        
        try:
            # 1. Load System (from a sanitized copy if its includes are indented)
            readable_top = GromacsPreprocessor.sanitize_topology(str(topology_file), output_dir)
            gmx = pmd.gromacs.GromacsTopologyFile(str(readable_top), xyz=str(coordinate_file))

            # 1b. Assign real GB radii (see docstring: GROMACS topologies have
            # none of their own) BEFORE saving, so every downstream .prmtop
            # carries correct, gb_model-matched radii rather than silently
            # defaulting to the crude element-only fallback later.
            try:
                from parmed.tools import changeRadii
                radii_set = GromacsPreprocessor._GB_MODEL_RADII_SET.get(gb_model, 'mbondi2')
                changeRadii(gmx, radii_set).execute()
                log.info(f"Applied '{radii_set}' GB radii set (matching gb_model='{gb_model}') "
                         f"to converted GROMACS structure.")
            except Exception as e:
                log.warning(f"Failed to apply GB radii set during Gromacs->Amber conversion: {e}. "
                            f"Resulting .prmtop files will have zero GB radii and downstream GB "
                            f"energies will NOT be comparable to Amber MMPBSA.py results.")

            # 2. Save Solvated Reference
            gmx.save(str(out_solvated), overwrite=True)
            
            # 3. Create Dry Complex
            dry_complex = gmx.copy(cls=pmd.Structure)
            dry_complex.strip(':SOL,NA,CL,TIP3,HOH,K,MG') # Standard water/ions
            dry_complex.save(str(out_dry), overwrite=True)
            
            # 4. Create Components (Receptor/Ligand)
            # Default assumption: Ligand is residue 'LIG' or detected
            # TODO: Better detection? For now, standard LIG.
            receptor = dry_complex.__copy__()
            receptor.strip(":LIG,UNL,UNK") 
            receptor.save(str(out_receptor), overwrite=True)
            
            ligand = dry_complex.__copy__()
            ligand.strip("!(:LIG,UNL,UNK)")
            ligand.save(str(out_ligand), overwrite=True)
            
            # 5. Fix Torsions (Charmm->Amber artifact)
            GromacsPreprocessor.fix_torsions([out_ligand, out_dry, out_receptor])
            
            log.info("Gromacs pre-processing complete.")
            
            return {
                'topology': str(out_dry),
                'solvated_topology': str(out_solvated),
                'receptor_topology': str(out_receptor),
                'ligand_topology': str(out_ligand)
            }
            
        except Exception as e:
            log.error(f"Gromacs conversion failed: {e}")
            raise e

    @staticmethod
    def sanitize_topology(topology_path, output_dir):
        """
        The path of a version of the GROMACS topology that ParmEd can read.

        ParmEd's preprocessor does not accept an indented #include. When the
        topology has one, a copy without the indentation is written to
        `output_dir`, with every relative #include that exists next to the
        original made absolute so that it is still found. The input file is
        never changed. Returns the original path when no change is needed.
        """
        topology_path = Path(topology_path)
        with open(topology_path, 'r') as f:
            lines = f.readlines()
        if not any(line.strip().startswith('#include') and line[:1].isspace() for line in lines):
            return topology_path

        new_lines = []
        for line in lines:
            stripped = line.strip()
            if stripped.startswith('#include'):
                quoted = stripped[len('#include'):].strip()
                name = quoted.strip('"<>')
                candidate = topology_path.parent / name
                if not Path(name).is_absolute() and candidate.exists():
                    stripped = f'#include "{candidate.resolve()}"'
                new_lines.append(stripped + '\n')
            else:
                new_lines.append(line)
        sanitized = Path(output_dir) / f"{topology_path.stem}_sanitized.top"
        with open(sanitized, 'w') as f:
            f.writelines(new_lines)
        log.info(f"Indented #include lines in {topology_path.name}: reading a sanitized copy, "
                 f"{sanitized}")
        return sanitized

    @staticmethod
    def fix_torsions(file_list):
        """Repairs invalid zero-periodicity dihedrals in prmtop files."""
        for f in file_list:
            f_path = Path(f)
            if not f_path.exists(): continue
            
            try:
                mol = pmd.load_file(str(f_path))
                changed = False
                dihedrals_to_remove = []
                processed_type_ids = set()
                
                # Iterate over all dihedrals
                for d in mol.dihedrals:
                    # Case 1: Single DihedralType
                    if isinstance(d.type, pmd.DihedralType):
                        if hasattr(d.type, 'per') and d.type.per <= 0:
                            dihedrals_to_remove.append(d)
                            changed = True
                    # Case 2: DihedralTypeList (Fourier Series)
                    elif isinstance(d.type, pmd.DihedralTypeList):
                        if id(d.type) in processed_type_ids:
                            continue
                        processed_type_ids.add(id(d.type))
                        
                        # Find bad terms in this list
                        bad_terms = [dt for dt in d.type if hasattr(dt, 'per') and dt.per <= 0]
                        if bad_terms:
                             for dt in bad_terms:
                                 d.type.remove(dt)
                             changed = True
                
                # Remove empty dihedrals (list became empty)
                for d in mol.dihedrals:
                    if isinstance(d.type, pmd.DihedralTypeList) and len(d.type) == 0:
                        dihedrals_to_remove.append(d)
                        changed = True

                # Remove Impropers with issues
                if hasattr(mol, 'impropers'):
                     for d in mol.impropers:
                        if isinstance(d.type, pmd.DihedralTypeList):
                            bad_terms = [dt for dt in d.type if hasattr(dt, 'per') and dt.per <= 0]
                            for dt in bad_terms:
                                d.type.remove(dt)
                                changed = True
                
                if dihedrals_to_remove:
                    log.info(f"Removing {len(dihedrals_to_remove)} invalid dihedrals from {f_path.name}")
                    for d in dihedrals_to_remove:
                        if d in mol.dihedrals:
                            mol.dihedrals.remove(d)
                
                if changed:
                    mol.prune_empty_terms()
                    mol.save(str(f_path), overwrite=True)
                    log.info(f"Fixed invalid torsions in {f_path.name}")
                    
            except Exception as e:
                log.warning(f"Could not check torsions for {f_path.name}: {e}")


class OpenMMSystemPreprocessor:
    """
    Converts a serialised OpenMM System (.xml) with a PDB file of the same
    atoms into an Amber .prmtop, so that it enters Native Mode like an Amber
    topology. The force-field parameters come from the System; the PDB gives
    the atom and residue names, the bonds and the coordinates. GB radii are not
    part of a System and are assigned for the GB model, as for GROMACS input.
    """

    # Forces that carry the force field and that ParmEd converts.
    CONVERTED = ('HarmonicBondForce', 'HarmonicAngleForce', 'PeriodicTorsionForce',
                 'NonbondedForce', 'CMAPTorsionForce')
    # Forces of the simulation setup that MM/GBSA does not use.
    DROPPED = ('CMMotionRemover', 'MonteCarloBarostat', 'MonteCarloAnisotropicBarostat',
               'MonteCarloMembraneBarostat', 'AndersenThermostat', 'GBSAOBCForce',
               'CustomGBForce')

    @staticmethod
    def find_pdb(xml_file):
        """The PDB that describes the atoms of `xml_file`: same name, .pdb suffix."""
        pdb = Path(xml_file).with_suffix('.pdb')
        if not pdb.exists():
            raise FileNotFoundError(
                f"An OpenMM System {Path(xml_file).name} needs a PDB file of the same atoms "
                f"with the same name, {pdb.name}, for the atom names, residues and bonds.")
        return pdb

    @staticmethod
    def convert_to_amber(xml_file, pdb_file, output_dir, gb_model='OBC2', radii_set='auto'):
        """Write <output_dir>/<stem>_from_openmm.prmtop and return its path.

        `radii_set` 'auto' assigns the set that `gb_model` was parameterised with
        (mbondi for HCT and GBn, mbondi2 for OBC1 and OBC2, mbondi3 for GBn2).
        """
        import openmm
        from openmm import app

        xml_file, pdb_file = Path(xml_file), Path(pdb_file)
        with open(xml_file) as fh:
            system = openmm.XmlSerializer.deserialize(fh.read())
        pdb = app.PDBFile(str(pdb_file))
        if pdb.topology.getNumAtoms() != system.getNumParticles():
            raise ValueError(f"{pdb_file.name} has {pdb.topology.getNumAtoms()} atoms but the "
                             f"System in {xml_file.name} has {system.getNumParticles()} particles.")

        for index in reversed(range(system.getNumForces())):
            name = type(system.getForce(index)).__name__
            if name in OpenMMSystemPreprocessor.DROPPED:
                system.removeForce(index)
            elif name not in OpenMMSystemPreprocessor.CONVERTED:
                raise ValueError(
                    f"{xml_file.name} contains a {name}, which cannot be converted to an Amber "
                    f"topology without changing the energy (for example the pair-specific "
                    f"Lennard-Jones terms of CHARMM). Give the original topology instead "
                    f"(.prmtop, .top or .psf).")

        # The bond graph must be the System's own. A PDB carries bonds from residue
        # templates and may join consecutive residues the System does not bond (a
        # protein and a nucleic acid in one chain, for example), and ParmEd checks
        # the 1-2, 1-3 and 1-4 exceptions of the NonbondedForce against the bonds.
        pairs = set()
        for force in system.getForces():
            if isinstance(force, openmm.HarmonicBondForce):
                for i in range(force.getNumBonds()):
                    a, b, *_ = force.getBondParameters(i)
                    pairs.add(tuple(sorted((a, b))))
        for i in range(system.getNumConstraints()):
            a, b, _ = system.getConstraintParameters(i)
            pairs.add(tuple(sorted((a, b))))
        topology = app.Topology()
        topology.setPeriodicBoxVectors(pdb.topology.getPeriodicBoxVectors())
        atoms = []
        for chain in pdb.topology.chains():
            new_chain = topology.addChain(chain.id)
            for residue in chain.residues():
                new_residue = topology.addResidue(residue.name, new_chain, residue.id,
                                                  residue.insertionCode)
                for atom in residue.atoms():
                    atoms.append(topology.addAtom(atom.name, atom.element, new_residue, atom.id))
        for a, b in sorted(pairs):
            topology.addBond(atoms[a], atoms[b])

        structure = pmd.openmm.load_topology(topology, system, xyz=pdb.positions)
        if structure.unknown_functional:
            raise ValueError(f"ParmEd could not express the System in {xml_file.name} as an Amber "
                             f"topology. Give the original topology instead (.prmtop, .top or .psf).")
        if not radii_set or radii_set == 'auto':
            radii_set = GromacsPreprocessor._GB_MODEL_RADII_SET.get(gb_model, 'mbondi2')
        from parmed.tools import changeRadii
        changeRadii(structure, radii_set).execute()
        log.info(f"Applied '{radii_set}' GB radii set to the converted OpenMM System "
                 f"(gb_model='{gb_model}').")

        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        out = output_dir / f"{xml_file.stem}_from_openmm.prmtop"
        structure.save(str(out), overwrite=True)
        return str(out)

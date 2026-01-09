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
    
    @staticmethod
    def convert_to_amber(topology_file, coordinate_file, output_dir=None):
        """
        Converts Gromacs topology to Dry Amber/Parm7 format.
        
        Returns:
            dict: Paths to generated files:
                - 'topology': complex_dry.prmtop (for simulation)
                - 'solvated': complex.prmtop (for trajectory mapping)
                - 'receptor': receptor.prmtop
                - 'ligand': ligand.prmtop
        """
        topology_file = Path(topology_file).resolve()
        coordinate_file = Path(coordinate_file).resolve()
        
        if output_dir is None:
            output_dir = topology_file.parent
        else:
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
            # 1. Load System
            # Pre-Sanitize to fix indented includes
            GromacsPreprocessor.sanitize_topology(str(topology_file))
            
            # Note: Requires ParmEd with indented include fix? 
            # We assume ParmEd is patched or input is clean.
            gmx = pmd.gromacs.GromacsTopologyFile(str(topology_file), xyz=str(coordinate_file))
            
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
    def sanitize_topology(topology_path):
        """
        Sanitizes Gromacs topology file to ensure compatibility with ParmEd.
        Specific fixes:
        1. Removes leading whitespace before #include handling (ParmEd parser issue).
        """
        try:
            with open(topology_path, 'r') as f:
                lines = f.readlines()
            
            changed = False
            new_lines = []
            for line in lines:
                # Fix indented includes
                if line.strip().startswith('#include') and line.startswith(' '):
                    new_lines.append(line.strip() + '\n')
                    changed = True
                else:
                    new_lines.append(line)
            
            if changed:
                log.info(f"Sanitizing topology: Fixed indented includes in {topology_path}")
                # Create backup
                backup = str(topology_path) + ".bak"
                if not os.path.exists(backup):
                    with open(backup, 'w') as f:
                        f.writelines(lines)
                
                # Overwrite with fixed content
                with open(topology_path, 'w') as f:
                    f.writelines(new_lines)
                    
        except Exception as e:
            log.warning(f"Failed to sanitize topology: {e}")

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

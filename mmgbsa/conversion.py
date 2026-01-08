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
    def fix_torsions(file_list):
        """Repairs invalid zero-periodicity dihedrals in prmtop files."""
        for f in file_list:
            f_path = Path(f)
            if not f_path.exists(): continue
            
            try:
                mol = pmd.load_file(str(f_path))
                changed = False
                for t in mol.dihedrals:
                    if t.type.per <= 0:
                        t.type.per = 1
                        changed = True
                
                if changed:
                    mol.save(str(f_path), overwrite=True)
                    log.info(f"Fixed invalid torsions in {f_path.name}")
            except Exception as e:
                log.warning(f"Could not check torsions for {f_path.name}: {e}")

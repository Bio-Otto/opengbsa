from pathlib import Path
from enum import Enum, auto

class EngineMode(Enum):
    AMBER = auto()      # .prmtop
    GROMACS = auto()    # .tpr
    CHARMM = auto()     # .psf
    OPENMM = auto()     # .xml
    GENERIC = auto()    # .pdb (Constructed)

class InputManager:
    """
    Handles detection of input file types and determines the simulation engine mode.
    """
    
    @staticmethod
    def detect_mode(topology_file: str) -> EngineMode:
        """
        Detects the engine mode based on the topology file extension.
        """
        path = Path(topology_file)
        ext = path.suffix.lower()
        
        if ext in ['.prmtop', '.parm7']:
            return EngineMode.AMBER
        elif ext in ['.tpr', '.top']:
            return EngineMode.GROMACS
        elif ext == '.psf':
            return EngineMode.CHARMM
        elif ext == '.xml':
            return EngineMode.OPENMM
        elif ext == '.pdb':
            return EngineMode.GENERIC
        else:
            raise ValueError(f"Unknown topology file extension: {ext}")

    @staticmethod
    def validate_inputs(config: dict, base_dir: Path = None) -> dict:
        """
        Validates and standardizes the configuration dictionary.
        base_dir: Optional Path to resolve relative headers against.
        Returns the processed config dictionary.
        """
        if 'input' not in config:
             raise ValueError("Configuration missing 'input' section.")

        inp = config['input']
        if 'topology' not in inp:
             raise ValueError("Configuration missing 'input.topology'.")
             
        # Resolve all paths
        for k in ['topology', 'trajectory', 'receptor_topology', 'ligand_topology', 'solvated_topology', 'ligand_mol', 'ligand_pdb']:
             if k in inp and inp[k]:
                  path_val = Path(inp[k])
                  if not path_val.is_absolute() and base_dir:
                       path_val = base_dir / path_val
                  inp[k] = str(path_val.resolve())
                  
        # Detect Mode
        mode = InputManager.detect_mode(inp['topology'])
        config['mode'] = mode
        
        return config

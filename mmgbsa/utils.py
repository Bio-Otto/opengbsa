"""
Utility Functions for MM/GBSA Analysis Package

This module contains utility functions used throughout the package.
"""

import os
import sys
import logging
import numpy as np
from pathlib import Path
from typing import Dict, Any, List, Optional, Tuple
import warnings

logger = logging.getLogger(__name__)

def setup_logging(level: str = "INFO", log_file: Optional[str] = None) -> None:
    """
    Setup logging configuration.
    
    Args:
        level: Logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL)
        log_file: Optional log file path
    """
    log_level = getattr(logging, level.upper(), logging.INFO)
    
    # Create formatter
    formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )
    
    # Setup root logger
    root_logger = logging.getLogger()
    root_logger.setLevel(log_level)
    
    # Clear existing handlers
    root_logger.handlers.clear()
    
    # Console handler
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(log_level)
    console_handler.setFormatter(formatter)
    root_logger.addHandler(console_handler)
    
    # File handler (if specified)
    if log_file:
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(log_level)
        file_handler.setFormatter(formatter)
        root_logger.addHandler(file_handler)
    
    logger.info(f"Logging setup complete. Level: {level}")

def create_output_directory(base_dir: str, analysis_name: str = "mmgbsa_analysis") -> str:
    """
    Create output directory with timestamp.
    
    Args:
        base_dir: Base directory path
        analysis_name: Analysis name prefix
        
    Returns:
        Created directory path
    """
    from datetime import datetime
    
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_dir = os.path.join(base_dir, f"{analysis_name}_{timestamp}")
    
    try:
        os.makedirs(output_dir, exist_ok=True)
        logger.info(f"Output directory created: {output_dir}")
        return output_dir
    except Exception as e:
        logger.error(f"Error creating output directory: {e}")
        raise

def validate_file_path(file_path: str, file_type: str = "file") -> bool:
    """
    Validate file path.
    
    Args:
        file_path: Path to validate
        file_type: Type of file (file, directory)
        
    Returns:
        True if valid, False otherwise
    """
    if not file_path:
        return False
    
    if file_type == "file":
        return os.path.isfile(file_path) and os.access(file_path, os.R_OK)
    elif file_type == "directory":
        return os.path.isdir(file_path) and os.access(file_path, os.R_OK)
    
    return False

def get_file_size_mb(file_path: str) -> float:
    """
    Get file size in megabytes.
    
    Args:
        file_path: Path to file
        
    Returns:
        File size in MB
    """
    try:
        size_bytes = os.path.getsize(file_path)
        return size_bytes / (1024 * 1024)
    except OSError:
        return 0.0

def format_time(seconds: float) -> str:
    """
    Format time in seconds to human readable string.
    
    Args:
        seconds: Time in seconds
        
    Returns:
        Formatted time string
    """
    if seconds < 60:
        return f"{seconds:.1f}s"
    elif seconds < 3600:
        minutes = seconds / 60
        return f"{minutes:.1f}m"
    else:
        hours = seconds / 3600
        return f"{hours:.1f}h"

def format_memory(bytes_value: int) -> str:
    """
    Format memory in bytes to human readable string.
    
    Args:
        bytes_value: Memory in bytes
        
    Returns:
        Formatted memory string
    """
    for unit in ['B', 'KB', 'MB', 'GB']:
        if bytes_value < 1024.0:
            return f"{bytes_value:.1f}{unit}"
        bytes_value /= 1024.0
    return f"{bytes_value:.1f}TB"

def calculate_statistics(data: List[float]) -> Dict[str, float]:
    """
    Calculate basic statistics for a list of values.
    
    Args:
        data: List of numerical values
        
    Returns:
        Dictionary with statistics
    """
    if not data:
        return {}
    
    data_array = np.array(data)
    
    return {
        'mean': float(np.mean(data_array)),
        'std': float(np.std(data_array)),
        'min': float(np.min(data_array)),
        'max': float(np.max(data_array)),
        'median': float(np.median(data_array)),
        'count': len(data_array)
    }

def safe_divide(numerator: float, denominator: float, default: float = 0.0) -> float:
    """
    Safely divide two numbers, returning default if denominator is zero.
    
    Args:
        numerator: Numerator
        denominator: Denominator
        default: Default value if division by zero
        
    Returns:
        Result of division or default value
    """
    try:
        return numerator / denominator if denominator != 0 else default
    except (TypeError, ZeroDivisionError):
        return default

def round_to_significant_figures(value: float, sig_figs: int = 3) -> float:
    """
    Round a number to a specified number of significant figures.
    
    Args:
        value: Number to round
        sig_figs: Number of significant figures
        
    Returns:
        Rounded number
    """
    if value == 0:
        return 0.0
    
    return round(value, sig_figs - int(np.floor(np.log10(abs(value)))) - 1)

def check_dependencies() -> Dict[str, bool]:
    """
    Check if required dependencies are available.
    
    Returns:
        Dictionary with dependency availability
    """
    dependencies = {
        'openmm': False,
        'mdtraj': False,
        'numpy': False,
        'pandas': False,
        'matplotlib': False,
        'seaborn': False,
        'yaml': False,
        'prolif': False,
        'rdkit': False,
        'MDAnalysis': False,
        'h5py': False,
        'openff-toolkit': False,
        'openmmforcefields': False
    }
    
    try:
        import openmm
        dependencies['openmm'] = True
    except ImportError:
        pass
    
    try:
        import mdtraj
        dependencies['mdtraj'] = True
    except ImportError:
        pass
    
    try:
        import numpy
        dependencies['numpy'] = True
    except ImportError:
        pass
    
    try:
        import pandas
        dependencies['pandas'] = True
    except ImportError:
        pass
    
    try:
        import matplotlib
        dependencies['matplotlib'] = True
    except ImportError:
        pass
    
    try:
        import seaborn
        dependencies['seaborn'] = True
    except ImportError:
        pass
    
    try:
        import yaml
        dependencies['yaml'] = True
    except ImportError:
        pass
    
    try:
        import prolif
        dependencies['prolif'] = True
    except ImportError:
        pass
    
    try:
        import rdkit
        dependencies['rdkit'] = True
    except ImportError:
        pass
    
    try:
        import MDAnalysis
        dependencies['MDAnalysis'] = True
    except ImportError:
        pass
    
    try:
        import h5py
        dependencies['h5py'] = True
    except ImportError:
        pass
    
    try:
        import openff
        dependencies['openff-toolkit'] = True
    except ImportError:
        pass
    
    try:
        import openmmforcefields
        dependencies['openmmforcefields'] = True
    except ImportError:
        pass
    
    return dependencies

def get_system_info() -> Dict[str, Any]:
    """
    Get system information.
    
    Returns:
        Dictionary with system information
    """
    import platform
    import psutil
    
    info = {
        'platform': platform.platform(),
        'python_version': platform.python_version(),
        'architecture': platform.architecture()[0],
        'processor': platform.processor(),
        'memory_total_gb': round(psutil.virtual_memory().total / (1024**3), 2),
        'cpu_count': psutil.cpu_count(),
        'cpu_count_logical': psutil.cpu_count(logical=True)
    }
    
    return info

def suppress_warnings():
    """Suppress common warnings."""
    warnings.filterwarnings('ignore', category=UserWarning)
    warnings.filterwarnings('ignore', category=DeprecationWarning)
    warnings.filterwarnings('ignore', category=FutureWarning)

def enable_warnings():
    """Enable all warnings."""
    warnings.resetwarnings()

def is_gpu_available() -> bool:
    """
    Check if GPU is available for OpenMM.
    
    Returns:
        True if GPU is available, False otherwise
    """
    try:
        import openmm
        platforms = openmm.Platform.getNumPlatforms()
        
        for i in range(platforms):
            platform = openmm.Platform.getPlatform(i)
            if platform.getName() in ['CUDA', 'OpenCL']:
                return True
        
        return False
    except ImportError:
        return False

def get_available_platforms() -> List[str]:
    """
    Get list of available OpenMM platforms.
    
    Returns:
        List of platform names
    """
    try:
        import openmm
        platforms = []
        
        for i in range(openmm.Platform.getNumPlatforms()):
            platform = openmm.Platform.getPlatform(i)
            platforms.append(platform.getName())
        
        return platforms
    except ImportError:
        return []

def estimate_memory_usage(n_atoms: int, n_frames: int) -> float:
    """
    Estimate memory usage for analysis.
    
    Args:
        n_atoms: Number of atoms
        n_frames: Number of frames
        
    Returns:
        Estimated memory usage in MB
    """
    # Rough estimation: 3 coordinates per atom per frame * 8 bytes per float
    memory_bytes = n_atoms * n_frames * 3 * 8
    return memory_bytes / (1024 * 1024)  # Convert to MB

def estimate_computation_time(n_atoms: int, n_frames: int, gb_model: str = "OBC2") -> float:
    """
    Estimate computation time for analysis.
    
    Args:
        n_atoms: Number of atoms
        n_frames: Number of frames
        gb_model: GB model being used
        
    Returns:
        Estimated time in seconds
    """
    # Rough estimation based on empirical data
    base_time_per_frame = 0.1  # seconds per frame for small systems
    
    # Scale with number of atoms (roughly quadratic)
    atom_factor = (n_atoms / 1000) ** 1.5
    
    # GB model factor
    gb_factors = {
        "OBC2": 1.0,
        "OBC1": 0.8,
        "HCT": 0.6,
        "GBn": 1.2,
        "GBn2": 1.4
    }
    gb_factor = gb_factors.get(gb_model, 1.0)
    
    estimated_time = base_time_per_frame * n_frames * atom_factor * gb_factor
    return estimated_time

def convert_mol2_to_sdf(mol2_path: str, sdf_path: str) -> bool:
    """
    Convert a Mol2 file to SDF using manual parsing and RDKit reconstruction.
    This bypasses RDKit's strict Mol2 parser and OpenFF's implicit protonation issues.
    
    Args:
        mol2_path: Path to input Mol2 file
        sdf_path: Path to output SDF file
        
    Returns:
        True if successful, False otherwise
    """
    import logging
    from rdkit import Chem
    
    logger = logging.getLogger(__name__)
    
    try:
        def parse_mol2_atom_line(line):
            parts = line.split()
            # atom_id, name, x, y, z, Type, ...
            # Mol2 format usually: atom_id atom_name x y z atom_type [subst_id [subst_name [charge [status_bit]]]]
            atom_id = int(parts[0])
            name = parts[1]
            x = float(parts[2])
            y = float(parts[3])
            z = float(parts[4])
            atom_type = parts[5]
            # Infer element from name or type
            # Simple heuristic: First one or two chars of type (e.g., C.3 -> C, N.ar -> N)
            symbol = atom_type.split('.')[0]
            # Handle numbers in symbol if any (rare in Mol2 types but possible)
            symbol = ''.join([c for c in symbol if c.isalpha()])
            return atom_id, symbol, (x, y, z)

        def parse_mol2_bond_line(line):
            parts = line.split()
            # bond_id, atom1, atom2, type
            atom1 = int(parts[1])
            atom2 = int(parts[2])
            bond_type_str = parts[3]
            
            if bond_type_str == '1':
                btype = Chem.BondType.SINGLE
            elif bond_type_str == '2':
                btype = Chem.BondType.DOUBLE
            elif bond_type_str == '3':
                btype = Chem.BondType.TRIPLE
            elif bond_type_str == 'ar':
                btype = Chem.BondType.AROMATIC
            elif bond_type_str == 'am': # Amide
                btype = Chem.BondType.SINGLE # Treat amide as single for now, RDKit handles resonance
            else:
                btype = Chem.BondType.SINGLE
                
            return atom1, atom2, btype

        with open(mol2_path, 'r') as f:
            lines = f.readlines()

        atoms = {}
        bonds = []
        
        section = None
        for line in lines:
            line = line.strip()
            if not line: continue
            if line.startswith('@<TRIPOS>'):
                section = line
                continue
                
            if section == '@<TRIPOS>ATOM':
                try:
                    atom_id, symbol, coords = parse_mol2_atom_line(line)
                    atoms[atom_id] = (symbol, coords)
                except (ValueError, IndexError):
                    continue
            elif section == '@<TRIPOS>BOND':
                try:
                    bonds.append(parse_mol2_bond_line(line))
                except (ValueError, IndexError):
                    continue

        if not atoms:
            logger.error(f"No atoms found in {mol2_path}")
            return False

        # Build RDKit Mol
        mw = Chem.RWMol()
        atom_mapping = {} # Mol2 ID -> RDKit Index
        
        sorted_ids = sorted(atoms.keys())
        for aid in sorted_ids:
            symbol, coords = atoms[aid]
            a = Chem.Atom(symbol)
            idx = mw.AddAtom(a)
            atom_mapping[aid] = idx

        # Add bonds
        for a1, a2, btype in bonds:
            if a1 in atom_mapping and a2 in atom_mapping:
                idx1 = atom_mapping[a1]
                idx2 = atom_mapping[a2]
                mw.AddBond(idx1, idx2, btype)

        mol = mw.GetMol()
        
        # Add Conformer
        conf = Chem.Conformer(mol.GetNumAtoms())
        for aid in sorted_ids:
            idx = atom_mapping[aid]
            coords = atoms[aid][1]
            conf.SetAtomPosition(idx, coords)
        mol.AddConformer(conf)
        
        try:
            Chem.SanitizeMol(mol)
        except Exception as e:
            logger.warning(f"Sanitization warning during Mol2 conversion: {e}")
            # Continue anyway, as we want to preserve atom count
        
        writer = Chem.SDWriter(sdf_path)
        # We might want to set Kekulize to False if aromaticity is tricky, 
        # but usually it's better to let RDKit try unless it fails.
        # If previous manual tests showed failure, we can set False.
        # But let's try default first, catching exception if needed.
        try:
            writer.write(mol)
        except Exception:
             writer.SetKekulize(False)
             writer.write(mol)
             
        writer.close()
        logger.info(f"Successfully converted {mol2_path} to {sdf_path} ({mol.GetNumAtoms()} atoms)")
        return True
        
    except Exception as e:
        logger.error(f"Error converting Mol2 to SDF: {e}")
        return False 
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
"""
MM/GBSA Analysis Package

A comprehensive Molecular Mechanics/Generalized Born Surface Area (MM/GBSA) 
analysis package with advanced features including entropy analysis, 
per-residue decomposition, and YAML-based configuration.

Author: H. Ibrahim Özdemir
License: MIT
Version: 0.0.2
"""

from ._version import __version__, __author__, __email__, __link__, __license__

# Core modules
from .runner import MMGBSARunner
# from .entropy import EntropyAnalyzer  # Removed: not defined
# from .decomposition import PerResidueDecomposer  # Removed: not defined
# from .visualization import MMGBSAVisualizer  # Removed: not defined

# Configuration
from .config import ConfigManager

__all__ = [
    "MMGBSARunner", 
    # "EntropyAnalyzer",  # Removed: not defined
    # "PerResidueDecomposer",  # Removed: not defined
    # "MMGBSAVisualizer",  # Removed: not defined
    "ConfigManager",
    "__version__",
    "__author__",
    "__email__",
    "__license__"
]

# Package metadata ve yardımcı fonksiyonlar burada kalabilir (isteğe bağlı) 
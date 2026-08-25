"""
Modular core MM/GBSA calculation components.

This package restructures the monolithic core.py into focused functional modules:
- platform: Platform configuration and management
- caching: System caching and retrieval
- topology: Structure building and atom selection
- parameterization: Ligand and protein parameterization
- analysis: Analysis execution and orchestration
  - FrameSelector: Multi-method trajectory sampling
  - EnergyCalculator: Frame-by-frame energy computation
  - DecompositionEngine: Per-residue energy analysis
  - ResultsAggregator: Statistics and validation
  - AnalysisEngine: High-level coordinator
- results: Result handling and reporting
  - ResultsExporter: Multi-format export (CSV/JSON)
  - ReportBuilder: HTML/text report generation
  - ResultsValidator: Quality control and validation
  - ResultsManager: Results orchestration
- calculator: Main GBSACalculator coordinator (facade)

Note: `StructureManager` and the low-level `GBSAForceManager`/system-building
code are not part of this package -- they live in `mmgbsa.mmgbsa_core`
(imported directly from there, not re-exported here).
"""

from .platform import PlatformManager
from .caching import CacheManager
from .topology import TopologyManager
from .parameterization import ParameterizationManager
from .analysis import (
    FrameSelector,
    EnergyCalculator,
    DecompositionEngine,
    ResultsAggregator,
    AnalysisEngine
)
from .results import (
    ResultsExporter,
    ReportBuilder,
    ResultsValidator,
    ResultsManager
)
from .calculator import GBSACalculator

__all__ = [
    # Platform management
    'PlatformManager',
    # Caching
    'CacheManager',
    # Topology
    'TopologyManager',
    # Parameterization
    'ParameterizationManager',
    # Analysis components
    'FrameSelector',
    'EnergyCalculator',
    'DecompositionEngine',
    'ResultsAggregator',
    'AnalysisEngine',
    # Results components
    'ResultsExporter',
    'ReportBuilder',
    'ResultsValidator',
    'ResultsManager',
    # Main facade
    'GBSACalculator',
]

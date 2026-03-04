"""
Facade GBSACalculator that composes new managers while delegating
to the legacy `mmgbsa_core.GBSACalculator` during incremental refactor.

This allows gradual extraction of functionality into managers while
keeping external API stable.
"""

from typing import Optional

from .platform import PlatformManager
from .caching import CacheManager
from .topology import TopologyManager
from .parameterization import ParameterizationManager
from .analysis import AnalysisEngine
from .results import ResultsManager


class GBSACalculator:
    """Facade that coordinates the refactored core components.

    Initially delegates heavy lifting to the legacy `mmgbsa_core.GBSACalculator`.
    """

    def __init__(self, temperature=300, verbose=1, gb_model='OBC2', salt_concentration=0.15, 
                 use_cache=True, parallel_processing=False, max_workers=None,
                 protein_forcefield='amber', charge_method='am1bcc', solute_dielectric=1.0,
                 solvent_dielectric=78.5, entropy_method='none', decomposition_method='full',
                 visualization_settings=None, platform=None, reporting_settings=None, 
                 sa_model='ACE', cache_dir=None, nonbonded_cutoff=None, **kwargs):
        """
        Initialize facade calculator with support for all legacy parameters.
        
        Parameters:
        -----------
        temperature : float
            Temperature in Kelvin
        verbose : int
            Verbosity level
        gb_model : str
            GB model (OBC1, OBC2, HCT, GBn, GBn2)
        salt_concentration : float
            Salt concentration in Molar
        use_cache : bool
            Enable system caching
        parallel_processing : bool
            Enable parallel frame processing
        max_workers : int, optional
            Max parallel workers
        protein_forcefield : str
            Protein forcefield
        charge_method : str
            Charge method for ligands (am1bcc, gasteiger)
        solute_dielectric : float
            Solute dielectric constant
        solvent_dielectric : float
            Solvent dielectric constant
        entropy_method : str
            Entropy calculation method
        decomposition_method : str
            Energy decomposition method
        visualization_settings : dict, optional
            Visualization configuration
        platform : str, optional
            Force specific platform
        reporting_settings : dict, optional
            Report generation settings
        sa_model : str
            Surface area model (ACE, LCPO)
        cache_dir : str, optional
            Cache directory path
        nonbonded_cutoff : float, optional
            Nonbonded interaction cutoff
        **kwargs : dict
            Additional legacy parameters
        """
        self.verbose = verbose
        self.temperature = temperature
        self.gb_model = gb_model
        self.salt_concentration = salt_concentration
        self.use_cache = use_cache
        self.parallel_processing = parallel_processing
        self.max_workers = max_workers
        self.protein_forcefield = protein_forcefield
        self.charge_method = charge_method
        self.solute_dielectric = solute_dielectric
        self.solvent_dielectric = solvent_dielectric
        self.entropy_method = entropy_method
        self.decomposition_method = decomposition_method
        self.visualization_settings = visualization_settings or {}
        self.platform = platform
        self.reporting_settings = reporting_settings or {}
        self.sa_model = sa_model
        self.cache_dir = cache_dir
        self.nonbonded_cutoff = nonbonded_cutoff
        self._extra_kwargs = kwargs

        # Legacy calculator (temporary) - lazy imported to avoid startup cost
        self._legacy = None

        # Managers
        self.platform_manager = PlatformManager(verbose=bool(verbose))
        self.cache = CacheManager(enabled=use_cache, verbose=bool(verbose))
        self.topology = TopologyManager(verbose=bool(verbose))
        self.parameterization = ParameterizationManager(verbose=bool(verbose))
        self.analysis = AnalysisEngine(calculator=self.legacy, verbose=bool(verbose))
        self.results = ResultsManager(verbose=bool(verbose))

    @property
    def legacy(self):
        """Lazy-load legacy calculator with all parameters."""
        if self._legacy is None:
            try:
                from .. import mmgbsa_core
                # Build kwargs for legacy calculator
                legacy_kwargs = {
                    'temperature': self.temperature,
                    'verbose': self.verbose,
                    'gb_model': self.gb_model,
                    'salt_concentration': self.salt_concentration,
                    'use_cache': self.use_cache,
                    'parallel_processing': self.parallel_processing,
                    'max_workers': self.max_workers,
                    'protein_forcefield': self.protein_forcefield,
                    'charge_method': self.charge_method,
                    'solute_dielectric': self.solute_dielectric,
                    'solvent_dielectric': self.solvent_dielectric,
                    'entropy_method': self.entropy_method,
                    'decomposition_method': self.decomposition_method,
                    'visualization_settings': self.visualization_settings,
                    'platform': self.platform,
                    'reporting_settings': self.reporting_settings,
                    'sa_model': self.sa_model,
                    'cache_dir': self.cache_dir,
                    'nonbonded_cutoff': self.nonbonded_cutoff,
                }
                # Add extra kwargs
                legacy_kwargs.update(self._extra_kwargs)
                
                self._legacy = mmgbsa_core.GBSACalculator(**legacy_kwargs)
            except Exception as e:
                if self.verbose:
                    print(f"Warning: Could not instantiate legacy calculator: {e}")
                self._legacy = None
        return self._legacy

    def set_platform_settings(self, platform_settings: Optional[dict]):
        """Configure platform preferences for analysis and decomposition."""
        self.platform_manager.set_platform_settings(platform_settings)
        if hasattr(self.analysis, 'calculator') and hasattr(self.analysis.calculator, 'set_platform_settings'):
            try:
                self.analysis.calculator.set_platform_settings(platform_settings)
            except Exception:
                pass
        # forward to legacy if present
        if self.legacy:
            try:
                self.legacy.set_platform_settings(platform_settings)
            except Exception:
                pass

    def run(self, *args, **kwargs):
        """Run analysis. Delegates to AnalysisEngine which may call legacy calculator."""
        # If analysis engine was constructed with legacy calculator, delegate
        return self.analysis.run(*args, **kwargs)

    def run_comprehensive(self, *args, **kwargs):
        return self.analysis.run_comprehensive(*args, **kwargs)

    def __getattr__(self, name):
        """
        Compatibility delegation: expose legacy calculator methods/attributes
        that are still used by decomposition/reporting modules.
        """
        legacy = self.legacy
        if legacy is not None and hasattr(legacy, name):
            return getattr(legacy, name)
        raise AttributeError(f"{self.__class__.__name__} has no attribute '{name}'")

    # Minimal compatibility methods
    def save_results(self, results, filename='results.json'):
        return self.results.save_results(results, filename=filename)

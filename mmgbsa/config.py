"""
Configuration Management Module for OpenGBSA

This module handles configuration file loading, validation, and management.
"""

import os
import yaml
from pathlib import Path
from typing import Dict, Any, Optional, List
import logging

logger = logging.getLogger(__name__)

class ConfigManager:
    """
    Configuration manager for MM/GBSA analysis.
    
    Handles loading, validation, and management of configuration files.
    """
    
    def __init__(self, config_path: Optional[str] = None):
        """
        Initialize configuration manager.
        
        Args:
            config_path: Path to configuration file
        """
        self.config_path = config_path
        self.config = {}
        self.validation_errors = []
        self.validation_warnings = []
        
        if config_path:
            self.load_config(config_path)
    
    def load_config(self, config_path: str) -> bool:
        """
        Load configuration from YAML file.
        
        Args:
            config_path: Path to configuration file
            
        Returns:
            True if successful, False otherwise
        """
        try:
            with open(config_path, 'r', encoding='utf-8') as f:
                self.config = yaml.safe_load(f)
            
            # Normalize configuration (handle legacy 'input' section)
            self._normalize_config()
            
            self.config_path = config_path
            logger.info(f"Configuration loaded from: {config_path}")
            return True
            
        except FileNotFoundError:
            logger.error(f"Configuration file not found: {config_path}")
            return False
        except yaml.YAMLError as e:
            logger.error(f"Error parsing YAML file: {e}")
            return False
        except Exception as e:
            logger.error(f"Error loading configuration: {e}")
            return False
    
    
    def _normalize_config(self):
        """
        Normalize configuration to standard format.
        Handles legacy 'input' section alias to 'input_files'.
        """
        if 'input_files' not in self.config and 'input' in self.config:
            logger.info("Normalizing legacy 'input' section to 'input_files'")
            inp = self.config['input']
            
            # Create input_files section
            self.config['input_files'] = {}
            
            # Map keys
            mapping = {
                'topology': 'complex_pdb',
                'trajectory': 'trajectory',
                'ligand_mol': 'ligand_mol',
                'ligand_pdb': 'ligand_pdb',
                'receptor_topology': 'receptor_topology',
                'ligand_topology': 'ligand_topology',
                'solvated_topology': 'solvated_topology'
            }
            
            for legacy_key, new_key in mapping.items():
                if legacy_key in inp:
                    self.config['input_files'][new_key] = inp[legacy_key]

    def validate_config(self) -> bool:
        """
        Validate configuration.
        
        Returns:
            True if valid, False otherwise
        """
        self.validation_errors = []
        self.validation_warnings = []
        
        # Check required sections
        required_sections = ['input_files', 'analysis_settings']
        for section in required_sections:
            if section not in self.config:
                self.validation_errors.append(f"Missing required section: {section}")
        
        if self.validation_errors:
            return False
        
        # Validate input files
        self._validate_input_files()
        
        # Validate analysis settings
        self._validate_analysis_settings()
        
        # Cross-field validation
        self._validate_cross_fields()
        
        return len(self.validation_errors) == 0
    
    def _validate_input_files(self):
        """Validate input files section."""
        input_files = self.config.get('input_files', {})
        
        required_files = ['complex_pdb', 'trajectory']
        for file_key in required_files:
            if file_key not in input_files:
                self.validation_errors.append(f"Missing required file: {file_key}")
                continue
            
            file_path = input_files[file_key]
            if not self._validate_file_path(file_path, file_key):
                self.validation_errors.append(f"Invalid file path: {file_path}")
    
    def _validate_analysis_settings(self):
        """Validate analysis settings section."""
        settings = self.config.get('analysis_settings', {})
        
        # Required parameters with validation rules
        required_params = {
            'temperature': (float, 0, 1000),
            'gb_model': (str, ['OBC2', 'OBC1', 'HCT', 'GBn', 'GBn2'], None),  # param_type, allowed_values, None
            'salt_concentration': (float, 0.0, 2.0),
            'max_frames': (int, 1, 10000)
        }
        
        for param, validation_rule in required_params.items():
            if param not in settings:
                # Inject defaults for common parameters instead of failing
                defaults = {
                    'temperature': 300.0,
                    'salt_concentration': 0.15,
                    'max_frames': 50
                }
                if param in defaults:
                    settings[param] = defaults[param]
                    self.config['analysis_settings'][param] = defaults[param] # Update actual config
                    continue
                else: 
                    self.validation_errors.append(f"Missing required parameter: {param}")
                    continue
            
            value = settings[param]
            
            if param == 'max_frames' and value is None:
                continue
                
            param_type, min_val, max_val = validation_rule
            
            if not self._validate_parameter(param, value, param_type, min_val, max_val):
                self.validation_errors.append(f"Invalid parameter value: {param} = {value}")
        
        # Validate binding_mode
        valid_binding_modes = ['standard', 'dimer_ligand', 'ppi']
        binding_mode = settings.get('binding_mode', 'standard')
        if binding_mode not in valid_binding_modes:
            self.validation_errors.append(
                f"Invalid binding_mode '{binding_mode}'. Must be one of: {valid_binding_modes}"
            )
        else:
            # Inject default so downstream code can always read it
            self.config['analysis_settings']['binding_mode'] = binding_mode
    
    def _validate_cross_fields(self):
        """Validate cross-field dependencies."""
        settings = self.config.get('analysis_settings', {})
        
        # Frame range validation
        frame_start = settings.get('frame_start')
        frame_end = settings.get('frame_end')
        if frame_start is not None and frame_end is not None:
            if frame_end <= frame_start:
                self.validation_errors.append("frame_end must be greater than frame_start")
        
        # Decomposition frames validation
        max_frames = settings.get('max_frames')
        decomp_frames = settings.get('decomp_frames')
        if decomp_frames is not None and max_frames is not None:
            if decomp_frames > max_frames:
                self.validation_errors.append("decomp_frames cannot be greater than max_frames")
        
        # Random seed validation
        frame_selection = settings.get('frame_selection')
        random_seed = settings.get('random_seed')
        if frame_selection == 'random' and random_seed is None:
            self.validation_errors.append("random_seed is required when frame_selection is 'random'")
    
    def _validate_file_path(self, file_path: str, file_key: str) -> bool:
        """Validate file path."""
        if not file_path:
            return False
        
        # Check if file exists
        if not os.path.exists(file_path):
            return False
        
        # Check if file is readable
        if not os.access(file_path, os.R_OK):
            return False
        
        # Check file size (max 1GB)
        file_size = os.path.getsize(file_path)
        if file_size > 1024 * 1024 * 1024:  # 1GB
            self.validation_warnings.append(f"Large file: {file_path} ({file_size / 1024**3:.1f}GB)")
        
        return True
    
    def _validate_parameter(self, param: str, value: Any, param_type: type, min_val: Any, max_val: Any) -> bool:
        """Validate parameter value."""
        # Type check
        if not isinstance(value, param_type):
            return False
        
        # Range check
        if isinstance(min_val, (int, float)) and isinstance(max_val, (int, float)):
            if value < min_val or value > max_val:
                return False
        elif isinstance(min_val, list):
            if value not in min_val:
                return False
        
        return True
    
    def get_config(self) -> Dict[str, Any]:
        """Get configuration dictionary."""
        return self.config.copy()
    
    def get_input_files(self) -> Dict[str, str]:
        """Get input files configuration."""
        return self.config.get('input_files', {}).copy()
    
    def get_analysis_settings(self) -> Dict[str, Any]:
        """Get analysis settings configuration."""
        return self.config.get('analysis_settings', {}).copy()
    
    def get_validation_errors(self) -> List[str]:
        """Get validation errors."""
        return self.validation_errors.copy()
    
    def get_validation_warnings(self) -> List[str]:
        """Get validation warnings."""
        return self.validation_warnings.copy()
    
    def create_default_config(self, output_path: str) -> bool:
        """
        Create default configuration file.
        
        Args:
            output_path: Path to save configuration file
            
        Returns:
            True if successful, False otherwise
        """
        default_config = {
            'input_files': {
                'ligand_mol': 'path/to/ligand.sdf',
                'complex_pdb': 'path/to/complex.pdb',
                'ligand_pdb': 'path/to/ligand.pdb',
                'trajectory': 'path/to/trajectory.xtc'
            },
            'analysis_settings': {
                'temperature': 300.0,
                'gb_model': 'OBC2',
                'salt_concentration': 0.15,
                'max_frames': 50,
                'frame_start': None,
                'frame_end': None,
                'frame_stride': None,
                'frame_selection': 'sequential',
                'random_seed': None,
                'run_entropy_analysis': False,
                'run_per_residue_decomposition': True,
                'decomp_frames': 10,
                'energy_decomposition': False,
                'use_cache': True,
                'parallel_processing': True,
                'use_gpu': False,
                'gpu_platform': None,
                'reimage_trajectory': True
            }
        }
        
        try:
            with open(output_path, 'w', encoding='utf-8') as f:
                yaml.dump(default_config, f, default_flow_style=False, indent=2)
            
            logger.info(f"Default configuration created: {output_path}")
            return True
            
        except Exception as e:
            logger.error(f"Error creating default configuration: {e}")
            return False
    
    def create_complete_config(self, output_path: str) -> bool:
        """
        Create complete configuration file with all options and explanations.
        
        Args:
            output_path: Path to save configuration file
            
        Returns:
            True if successful, False otherwise
        """
        complete_config_yaml = """# OpenGBSA Complete Configuration File
# This file contains all available settings with explanations and alternative options.

input_files:
  # Path to the ligand file (SDF format is highly recommended over PDB for preserving bond orders and formal charges)
  ligand_mol: 'test/ligand.sdf'
  
  # Path to the complex PDB file (receptor + ligand)
  complex_pdb: 'test/complex.pdb'
  
  # Path to the ligand-only PDB file
  ligand_pdb: 'test/ligand.pdb'
  
  # Path to the MD trajectory file (XTC, DCD, TRR, etc.)
  trajectory: 'test/complex.xtc'
  
  # Optional: Explicit Gromacs/Amber topology files (leave commented if using PDB)
  # receptor_topology: 'path/to/receptor.top'
  # ligand_topology: 'path/to/ligand.top'
  # solvated_topology: 'path/to/complex_solvated.top'

output_settings:
  # Main directory where all analysis results will be saved
  output_directory: 'mmgbsa_results'
  
  # Name of this specific analysis run (used for subdirectories)
  analysis_name: 'sample_analysis'
  
  # Output formats to save results in
  output_formats: ['csv', 'txt', 'yaml']
  
  # Enable saving plots automatically
  save_plots: true
  
  # Formats for saving plots
  plot_formats: ['png', 'pdf']
  
  # Save intermediate generated files (useful for debugging, takes more space)
  save_intermediate: false
  
  # Save aligned/processed trajectories
  save_trajectories: false
  
  # Save detailed log files
  save_logs: true
  
  # Compress large output files automatically
  compress_output: false

analysis_settings:
  # === Core Execution Settings ===
  
  # Analytical binding mode. Options:
  # - 'standard': Standard Receptor-Ligand binding.
  # - 'dimer_ligand': Dimer-Ligand binding (Ligand is extracted from dimer interface).
  # - 'ppi': Protein-Protein Interaction (Treats one protein chain as the ligand).
  binding_mode: 'standard'
  
  # Temperature in Kelvin for calculations
  temperature: 310.0
  
  # Generalized Born (GB) model for solvation free energy. Options:
  # - 'OBC2' (Default, recommended for proteins)
  # - 'OBC1' 
  # - 'HCT'
  # - 'GBn'
  # - 'GBn2'
  gb_model: 'OBC2'
  
  # Salt concentration in Molar (M) for the implicit solvent model
  salt_concentration: 0.15

  # === Frame Selection Strategy ===
  
  # Maximum number of frames to process
  max_frames: 100
  
  # Starting frame index (0-indexed). If null, starts from the beginning.
  frame_start: null
  
  # Ending frame index. If null, goes to the end of the trajectory.
  frame_end: null
  
  # Stride (step size) between frames. Use 1 for every frame, 10 for every 10th frame.
  frame_stride: null
  
  # How to select frames if the trajectory has more than max_frames. Options:
  # - 'sequential': Selects frames evenly spaced across the trajectory.
  # - 'random': Selects frames randomly (requires random_seed).
  frame_selection: 'sequential'
  
  # Random seed for reproducible random frame selection. Used only if frame_selection is 'random'.
  random_seed: 42

  # Re-image molecules across periodic boundaries (PBC re-wrapping) before
  # computing any energy. A raw MD trajectory has no guarantee a molecule
  # stays whole/centered from frame to frame -- if the protein or ligand
  # drifts and wraps to the opposite side of the periodic box mid-trajectory,
  # downstream vdW/electrostatic energies are silently corrupted even though
  # the true physics hasn't changed. True/False. Default true. Only disable
  # if the input trajectory is already known to be correctly imaged (e.g.
  # already processed with cpptraj `autoimage` or an equivalent tool) --
  # re-imaging an already-imaged trajectory is a harmless no-op but still
  # costs time, so this lets you skip it in that specific case.
  reimage_trajectory: true

  # === Advanced Analysis Options ===
  
  # Calculate entropy (Interaction Entropy method). True/False.
  # Note: Entropy calculations can be computationally expensive.
  run_entropy_analysis: true
  
  # Perform per-residue energy decomposition. True/False.
  run_per_residue_decomposition: true
  
  # Number of frames to use for per-residue decomposition. Must be <= max_frames.
  decomp_frames: 20
  
  # Enable pairwise energy decomposition (Residue-Residue interaction matrix). True/False.
  energy_decomposition: true
  
  # Save frame-by-frame decomposition energies to CSV
  save_frame_by_frame_csv: true
  
  # Filename for the frame-by-frame CSV
  frame_by_frame_csv_name: "frame_by_frame_decomposition"
  
  # Include overall residue averages in the output summaries
  include_residue_summary: true
  
  # Components to output in the frame-by-frame file
  frame_output_components: ['vdw', 'electrostatic', 'solvation', 'total']
  
  # Format for detailed decomposition outputs (csv, json, hdf5)
  frame_output_format: 'csv'

  # === Hardware & Performance ===
  
  # Cache intermediate topologies and calculations to speed up re-runs. True/False.
  use_cache: true
  
  # Use Python multiprocessing to calculate frames in parallel (CPU only). True/False.
  parallel_processing: true
  
  # Maximum number of parallel workers (leave null for auto-detect based on CPU cores)
  max_workers: null

forcefield_settings:
  # Protein forcefield to use for parameterization
  protein_forcefield: 'amber14-all.xml'
  
  # Specific variant of the protein forcefield (optional)
  protein_variant: null
  
  # Ligand parameterization forcefield
  ligand_forcefield: 'openff-2.1.0.offxml'
  
  # Water model to use
  water_model: null
  
  # Ion parameters to use
  ion_forcefield: null

advanced_settings:
  # Tolerance for OpenMM energy minimization
  minimization_tolerance: 1.0e-06
  
  # Maximum steps for energy minimization before analysis
  max_minimization_steps: 10000
  
  # Normal Mode Analysis quality threshold
  nma_quality_threshold: 'Good'
  
  # Hot spot interaction energy threshold (kcal/mol). Residues below this are flagged as hot spots.
  hot_spot_threshold: -1.0
  
  # Bootstrap confidence interval (e.g. 0.95 for 95%)
  bootstrap_confidence: 0.95
  
  # Number of bootstrap samples for error estimation
  bootstrap_samples: 1000

platform_settings:
  # Preferred OpenMM platform. Options: 'CUDA', 'OpenCL', 'CPU', 'Reference'
  preferred_platform: 'CUDA'
  
  # CUDA device index to use if platform is CUDA
  cuda_device: 0
  
  # Precision for CUDA calculations. Options: 'single', 'mixed', 'double'
  cuda_precision: 'mixed'

reporting_settings:
  # Generate a comprehensive final report
  generate_final_report: true
  
  # Automatically generate result plots
  include_plots: true
  
  # Include summary of configuration in the final report
  include_config_summary: true

reproducibility_settings:
  # Save a copy of the configuration next to the results
  save_configuration: true
  
  # Save environment info to track Python/package versions
  save_environment: true
"""
        
        try:
            with open(output_path, 'w', encoding='utf-8') as f:
                f.write(complete_config_yaml)
            
            logger.info(f"Complete configuration created: {output_path}")
            return True
            
        except Exception as e:
            logger.error(f"Error creating complete configuration: {e}")
            return False 
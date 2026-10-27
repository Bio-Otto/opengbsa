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
        Handles legacy 'input' section alias to 'input_files', and turns a
        section left empty in YAML (which loads as null) into an empty section.
        """
        for section, value in list(self.config.items()):
            if value is None and section.endswith(('_files', '_settings', 'params')):
                self.config[section] = {}
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
                    'gb_model': 'OBC2',
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
            
            if param_type is float and isinstance(value, int) and not isinstance(value, bool):
                value = float(value)   # 300 means 300.0
                self.config['analysis_settings'][param] = value
            if not self._validate_parameter(param, value, param_type, min_val, max_val):
                self.validation_errors.append(f"Invalid parameter value: {param} = {value}")

        # Options that select a code path: an unknown value would otherwise be ignored
        # or fail only later in the run
        pff = (self.config.get('forcefield_settings') or {}).get('protein_forcefield')
        allowed_pff = ['amber', 'amber14', 'amber14-all', 'amber14-all.xml', 'charmm',
                       'charmm_gromacs', 'auto']
        if pff is not None and str(pff).strip().lower() not in allowed_pff:
            self.validation_errors.append(f"Invalid protein_forcefield '{pff}'. Must be one of: "
                                          f"['amber', 'amber14', 'charmm', 'charmm_gromacs']")
        for param, allowed in (('frame_selection', ['sequential', 'equidistant', 'random']),
                               ('entropy_method', ['none', 'interaction', 'c2', 'quasiharmonic',
                                                   'normal_mode', 'nmode']),
                               ('decomposition_method', ['full', 'fast'])):
            value = settings.get(param)
            if value is not None and value not in allowed:
                self.validation_errors.append(f"Invalid {param} '{value}'. Must be one of: {allowed}")
        
        # Validate gb_radii_set ('auto': an Amber .prmtop keeps its own radii)
        valid_radii_sets = ['auto', 'mbondi', 'mbondi2', 'mbondi3', 'bondi']
        radii_set = settings.get('gb_radii_set', 'auto')
        if radii_set not in valid_radii_sets:
            self.validation_errors.append(
                f"Invalid gb_radii_set '{radii_set}'. Must be one of: {valid_radii_sets}"
            )

        # Born-radius cutoff (A, Amber's rgbmax) and the MMPBSA.py-compatible preset
        rgbmax = settings.get('rgbmax', None)
        if rgbmax is not None and (not isinstance(rgbmax, (int, float)) or isinstance(rgbmax, bool) or rgbmax <= 0):
            self.validation_errors.append(f"Invalid rgbmax '{rgbmax}': a positive distance in A, or null for no cutoff")
        if not isinstance(settings.get('amber_compatible', False), bool):
            self.validation_errors.append("amber_compatible must be true or false")
        lcpo_typing = settings.get('lcpo_typing', None)
        if lcpo_typing is not None and lcpo_typing not in ('element', 'amber'):
            self.validation_errors.append(f"Invalid lcpo_typing '{lcpo_typing}'. Must be one of: ['element', 'amber']")
        sa_model = settings.get('sa_model', None)
        if sa_model is not None and sa_model not in ('ACE', 'LCPO'):
            self.validation_errors.append(f"Invalid sa_model '{sa_model}'. Must be one of: ['ACE', 'LCPO']")

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
        
        # Decomposition frames validation (only when a decomposition runs)
        max_frames = settings.get('max_frames')
        decomp_frames = settings.get('decomp_frames')
        decomposition_on = (settings.get('run_per_residue_decomposition') or
                            settings.get('per_residue_decomposition') or
                            (self.config.get('advanced_settings') or {}).get('run_per_residue_decomposition'))
        if decomposition_on and decomp_frames is not None and max_frames is not None:
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
        Create a short, commented configuration file with the settings most
        runs need (`opengbsa --create-config`).

        Args:
            output_path: Path to save configuration file

        Returns:
            True if successful, False otherwise
        """
        default_config_yaml = """# OpenGBSA configuration (opengbsa --create-config).
# Run it with:  opengbsa mmgbsa_config.yaml
# Every option is listed with --create-complete-config and explained at
# https://opengbsa.readthedocs.io/

input_files:
  # Native Mode: the topology used for the MD simulation (.prmtop, .top, .tpr,
  # .psf or a serialised OpenMM system .xml with a .pdb of the same name).
  # Coordinate Mode: a .pdb of the complex; OpenGBSA then builds the topology.
  complex_pdb: path/to/complex.prmtop
  trajectory: path/to/trajectory.xtc
  # Coordinate Mode only: the small-molecule ligand with bond orders and charges.
  # ligand_mol: path/to/ligand.sdf
  # Native Mode, optional: separate receptor and ligand topologies.
  # receptor_topology: path/to/receptor.prmtop
  # ligand_topology: path/to/ligand.prmtop
  # Topology of the atoms in the trajectory, if it still contains solvent.
  # solvated_topology: path/to/complex_solvated.gro

output_settings:
  output_directory: mmgbsa_results

analysis_settings:
  gb_model: OBC2              # OBC2, OBC1, HCT, GBn or GBn2
  salt_concentration: 0.15    # M
  temperature: 300.0          # K
  max_frames: 50
  frame_start: null           # first frame (0-based); null = start
  frame_end: null             # end frame (exclusive); null = last
  frame_stride: null          # every Nth frame; null = 1
  run_per_residue_decomposition: false
  decomp_frames: 10
  run_entropy_analysis: false
  entropy_method: interaction # interaction or c2
  amber_compatible: false     # true = reproduce Amber MMPBSA.py's conventions
"""

        try:
            with open(output_path, 'w', encoding='utf-8') as f:
                f.write(default_config_yaml)

            logger.info(f"Default configuration created: {output_path}")
            return True

        except Exception as e:
            logger.error(f"Error creating default configuration: {e}")
            return False

    def create_complete_config(self, output_path: str) -> bool:
        """
        Create a configuration file with every option the analysis reads, its
        default and its meaning (`opengbsa --create-complete-config`).

        Args:
            output_path: Path to save configuration file

        Returns:
            True if successful, False otherwise
        """
        complete_config_yaml = """# OpenGBSA complete configuration (opengbsa --create-complete-config).
# Every option the analysis reads, set to its default. Delete what you do not
# need: an option that is left out takes the value shown here.
# Reference: https://opengbsa.readthedocs.io/

# ---------------------------------------------------------------- input files
input_files:
  # Native Mode (the topology of the MD simulation is read as is):
  #   Amber .prmtop, GROMACS .top (with a .gro of the same folder) or .tpr,
  #   CHARMM/NAMD .psf (needs forcefield_settings.charmm_params), or a
  #   serialised OpenMM System .xml with a .pdb of the same name next to it.
  # Coordinate Mode (OpenGBSA builds the topology from a force field):
  #   .pdb or .gro of the complex, plus ligand_mol for a small-molecule ligand.
  complex_pdb: path/to/complex.prmtop

  # MD trajectory (.xtc, .trr, .dcd, .nc, .trj, multi-model .pdb, ...).
  trajectory: path/to/trajectory.xtc

  # Topology of the atoms in the trajectory when it still contains water and
  # ions (for example the solvated .gro or .prmtop). Not needed for a dry
  # trajectory.
  solvated_topology: null

  # Native Mode, Amber only: receptor and ligand topologies written with the
  # complex. Their atom counts define the receptor and the ligand.
  receptor_topology: null
  ligand_topology: null

  # Coordinate Mode: the small-molecule ligand (.sdf, or .mol2), with bond
  # orders and partial charges. Parameterised with OpenFF.
  ligand_mol: null

  # Optional ligand-only .pdb; used only to check the ligand residue names.
  ligand_pdb: null

  # Residue name shown in the ligand box of the interactive 3D viewer.
  ligand_resname: null

# ---------------------------------------------------------------- output
output_settings:
  # Each run writes <output_directory>/analysis_<date>_<time>/.
  # opengbsa --output-dir DIR overrides this.
  output_directory: mmgbsa_results

# ---------------------------------------------------------------- analysis
analysis_settings:
  # --- receptor and ligand ---
  # standard: receptor = everything but the ligand, ligand found automatically
  #   (receptor_topology atom count, else residue LIG/UNL/UNK/DRUG/MOL, else the
  #   smallest non-solvent residue). ppi: receptor = chainid 0, ligand = chainid 1.
  # dimer_ligand: same as standard.
  binding_mode: standard
  # MDTraj selections; they override binding_mode, e.g. "chainid 1".
  receptor_selection: null
  ligand_selection: null

  # --- solvation model ---
  gb_model: OBC2               # OBC2, OBC1, HCT, GBn, GBn2
  gb_radii_set: auto           # auto, mbondi, mbondi2, mbondi3, bondi
  salt_concentration: 0.15     # M, Debye-Hueckel screening
  temperature: 300.0           # K
  solute_dielectric: 1.0
  solvent_dielectric: null     # null = 78.5 (MMPBSA.py's value with amber_compatible)
  sa_model: LCPO               # non-polar term: LCPO or ACE
  lcpo_typing: null            # element (default) or amber
  rgbmax: null                 # Born-radius cutoff in A; null = every atom pair
  nonbonded_cutoff: null       # A; null (or >= 999) = no cutoff
  # true = Amber MMPBSA.py's conventions: rgbmax 25 A, MMPBSA.py's solvent
  # dielectric and LCPO with Amber's atom typing.
  amber_compatible: false

  # --- frames ---
  frame_start: null            # 0-based; null = first frame
  frame_end: null              # exclusive; null = last frame
  frame_stride: null           # null = 1
  max_frames: 50               # null = every selected frame
  # sequential: the first max_frames frames of the range after striding;
  # equidistant: max_frames frames spread evenly over the range;
  # random: max_frames random frames (random_seed).
  frame_selection: sequential
  random_seed: 42
  reimage_trajectory: true     # make molecules whole across periodic boundaries

  # --- per-residue decomposition ---
  run_per_residue_decomposition: false
  decomp_frames: 10            # frames used; must not exceed max_frames
  report_raw_energies: false   # also write complex/receptor/ligand energies
  plot_top_residues: 10        # residues shown in the hot-spot plots
  # A second, faster decomposition inside the main calculation (fast or full),
  # mainly for development.
  energy_decomposition: false
  decomposition_method: full

  # --- entropy (-TdS) ---
  # run_entropy_analysis switches the estimate on; entropy_method chooses it:
  # interaction (interaction entropy) or c2 (second-order cumulant), both from
  # the fluctuation of the gas-phase interaction energy and reported with its
  # standard deviation (IE unreliable above 3.6, C2 above 6.0 kcal/mol);
  # quasiharmonic: ligand quasi-harmonic entropy (qha_analyze_complex: true
  # adds receptor and complex to give a binding entropy); normal_mode:
  # normal-mode entropy of the ligand alone (needs ligand_mol), not added to dG.
  run_entropy_analysis: false
  entropy_method: interaction
  qha_analyze_complex: false
  nmode_cutoff_frequency: 10.0 # normal_mode: modes below this (cm-1) are left out

  # --- execution ---
  max_workers: null            # CPU worker processes; null = min(cores - 1, 4)
  print_interval: 10           # progress line every N frames
  use_cache: true              # reuse parameterised systems within a run

# ---------------------------------------------------------------- force field
forcefield_settings:
  # Coordinate Mode only: amber (ff14SB, default), amber14 (amber14-all),
  # charmm or charmm_gromacs (CHARMM36).
  protein_forcefield: amber
  # CHARMM/NAMD .psf input: parameter files (.prm, .str, .rtf), one path or a list.
  charmm_params: null
  # CHARMM/NAMD .psf input: coordinates (.pdb, .coor, .crd) if no
  # solvated_topology is given.
  charmm_coordinates: null

# ---------------------------------------------------------------- hardware
platform_settings:
  # null = automatic (CUDA, then OpenCL, then CPU). Or CUDA, OpenCL, CPU,
  # Reference. Select a GPU with the environment variable CUDA_VISIBLE_DEVICES.
  preferred_platform: null
  # Platform for the per-residue decomposition; null = same as above.
  decomposition_platform: null
  # true = stop with an error instead of falling back to another platform.
  strict_platform: false

advanced_settings:
  # Worker processes for a per-residue decomposition on the CPU; -1 = all cores.
  n_jobs: -1

# ---------------------------------------------------------------- HTML report
reporting_settings:
  include_all_plots: true      # energy time series and distributions
  export_format: png           # format of the plot export button
  export_width: 4200
  export_height: 2800
"""

        try:
            with open(output_path, 'w', encoding='utf-8') as f:
                f.write(complete_config_yaml)
            
            logger.info(f"Complete configuration created: {output_path}")
            return True
            
        except Exception as e:
            logger.error(f"Error creating complete configuration: {e}")
            return False 
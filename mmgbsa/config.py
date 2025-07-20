"""
Configuration Management Module for MM/GBSA Analysis Package

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
        
        required_files = ['ligand_mol', 'complex_pdb', 'ligand_pdb', 'trajectory']
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
                self.validation_errors.append(f"Missing required parameter: {param}")
                continue
            
            value = settings[param]
            param_type, min_val, max_val = validation_rule
            
            if not self._validate_parameter(param, value, param_type, min_val, max_val):
                self.validation_errors.append(f"Invalid parameter value: {param} = {value}")
    
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
        max_frames = settings.get('max_frames', 0)
        decomp_frames = settings.get('decomp_frames', 0)
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
                'gpu_platform': None
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
        Create complete configuration file with all options.
        
        Args:
            output_path: Path to save configuration file
            
        Returns:
            True if successful, False otherwise
        """
        complete_config = {
            'input_files': {
                'ligand_mol': 'path/to/ligand.sdf',
                'complex_pdb': 'path/to/complex.pdb',
                'ligand_pdb': 'path/to/ligand.pdb',
                'trajectory': 'path/to/trajectory.xtc'
            },
            'analysis_settings': {
                # Core parameters
                'temperature': 300.0,
                'gb_model': 'OBC2',
                'salt_concentration': 0.15,
                
                # Frame selection
                'max_frames': 50,
                'frame_start': None,
                'frame_end': None,
                'frame_stride': None,
                'frame_selection': 'sequential',
                'random_seed': None,
                
                # Analysis options
                'run_entropy_analysis': True,
                'run_per_residue_decomposition': True,
                'decomp_frames': 10,
                'energy_decomposition': True,
                
                # Performance
                'use_cache': True,
                'parallel_processing': True,
                'use_gpu': False,
                'gpu_platform': None
            }
        }
        
        try:
            with open(output_path, 'w', encoding='utf-8') as f:
                yaml.dump(complete_config, f, default_flow_style=False, indent=2)
            
            logger.info(f"Complete configuration created: {output_path}")
            return True
            
        except Exception as e:
            logger.error(f"Error creating complete configuration: {e}")
            return False 
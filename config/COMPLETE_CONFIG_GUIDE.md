# Complete MM/GBSA Configuration Guide

## 🎯 Overview

The complete configuration system provides full control over all aspects of MM/GBSA analysis with maximum reproducibility and flexibility.

## 📁 Files

- **`complete_mmgbsa_config.yaml`**: Complete configuration with all parameters and alternatives
- **`complete_mmgbsa_runner.py`**: Runner that handles the complete configuration
- **`test_complete_config.py`**: Test script to validate configuration

## 🚀 Quick Start

```bash
# 1. Create complete configuration
python complete_mmgbsa_runner.py --create-config

# 2. Edit configuration file
nano complete_mmgbsa_config.yaml

# 3. Run analysis
python complete_mmgbsa_runner.py complete_mmgbsa_config.yaml

# 4. Test configuration
python test_complete_config.py
```

## 📋 Configuration Sections

### 1. Input Files (`input_files`)
```yaml
input_files:
  ligand_mol: "test/ligand.sdf"      # Required: Ligand structure
  complex_pdb: "test/complex.pdb"    # Required: Protein-ligand complex
  ligand_pdb: "test/ligand.pdb"      # Required: Isolated ligand
  trajectory: "test/complex.xtc"     # Required: MD trajectory
  reference_structure: null          # Optional: For alignment
  additional_ligands: null           # Optional: For comparison
```

**File Format Alternatives:**
- Ligand: `.sdf`, `.mol2`, `.pdb`, `.xyz`
- Complex: `.pdb`, `.pdbqt`, `.gro`
- Trajectory: `.xtc`, `.dcd`, `.trr`, `.nc`, `.h5`

### 2. Output Settings (`output_settings`)
```yaml
output_settings:
  output_directory: "mmgbsa_results"     # Main output directory
  analysis_name: "complete_analysis"      # Analysis subdirectory name
  output_formats: ["csv", "txt", "yaml"] # Output file formats
  save_plots: true                       # Generate plots
  plot_formats: ["png", "pdf"]           # Plot file formats
  save_intermediate: false               # Save intermediate files
  save_trajectories: false               # Save processed trajectories
  save_logs: true                        # Save detailed logs
  compress_output: false                 # Compress large files
```

### 3. Analysis Settings (`analysis_settings`)
```yaml
analysis_settings:
  # Core parameters
  temperature: 300.0                    # Temperature (K)
  gb_model: "OBC2"                      # GB model: OBC2, OBC1, HCT, GBn, GBn2
  salt_concentration: 0.15              # Salt concentration (M)
  
  # Frame selection
  max_frames: 50                        # Maximum frames to analyze
  frame_start: null                     # Start frame (null = beginning)
  frame_end: null                       # End frame (null = end)
  frame_stride: null                    # Frame stride (null = all frames)
  frame_selection: "sequential"         # sequential, equidistant, random
  random_seed: 42                       # For random selection
  
  # Analysis options
  verbose: 1                            # Verbosity level (0-3)
  use_cache: true                       # Enable caching
  parallel_processing: false            # Enable parallel processing
  max_workers: null                     # Number of workers (null = auto)
  energy_decomposition: false           # Enable energy decomposition
  run_entropy_analysis: false           # Run entropy analysis
  run_per_residue_decomposition: false  # Run per-residue analysis
  decomp_frames: 10                     # Frames for decomposition
```

### 4. Forcefield Settings (`forcefield_settings`)
```yaml
forcefield_settings:
  # Protein forcefield
  protein_forcefield: "amber14-all.xml"     # amber14-all.xml, amber99sb.xml, charmm36.xml
  protein_variant: null                      # ff14SB, ff99SB, ff03, null
  
  # Ligand forcefield (OpenFF)
  ligand_forcefield: "openff-2.1.0.offxml"  # openff-2.1.0.offxml, openff-2.0.0.offxml
  ligand_variant: null                       # sage, parsley, smirnoff99Frosst, null
  
  # Water model
  water_model: null                          # tip3p.xml, tip4pew.xml, spce.xml, null
  
  # Ion forcefield
  ion_forcefield: null                       # amber14-all.xml, charmm36.xml, null
```

### 5. Advanced Settings (`advanced_settings`)
```yaml
advanced_settings:
  # Minimization
  minimization_tolerance: 1e-6              # Minimization tolerance
  max_minimization_steps: 10000             # Maximum minimization steps
  
  # Normal Mode Analysis
  nma_quality_threshold: "Good"             # Poor, Fair, Good, Excellent
  nma_tweak_ratio: 1e-12                    # NMA tweak energy ratio
  
  # Hot Spot Analysis
  hot_spot_threshold: -1.0                  # Hot spot energy threshold (kcal/mol)
  
  # Statistical Analysis
  bootstrap_confidence: 0.95                # Bootstrap confidence level
  bootstrap_samples: 1000                   # Number of bootstrap samples
  convergence_window: 10                    # Convergence window size
  
  # Caching
  cache_directory: ".cache"                 # Cache directory
  max_cache_size: 1000                      # Maximum cache size (MB)
  cache_expiration: 30                      # Cache expiration (days)
```

### 6. Platform Settings (`platform_settings`)
```yaml
platform_settings:
  preferred_platform: "CUDA"                # CUDA, OpenCL, CPU, Reference
  cuda_device: 0                           # CUDA device number
  cuda_precision: "double"                 # single, double, mixed
  deterministic_forces: true               # Enable deterministic forces
  platform_properties: {}                  # Additional platform properties
```

### 7. Validation Settings (`validation_settings`)
```yaml
validation_settings:
  validate_inputs: true                    # Validate input files
  validate_results: true                   # Validate results
  max_std_dev: 10.0                        # Maximum allowed std dev (kcal/mol)
  min_successful_frames: 10                # Minimum successful frames
  check_convergence: true                  # Check for convergence
  generate_validation_report: true         # Generate validation report
  validation_tolerance: 1.0                # Validation tolerance
```

### 8. Reporting Settings (`reporting_settings`)
```yaml
reporting_settings:
  generate_final_report: true              # Generate comprehensive report
  include_plots: true                      # Include plots in report
  include_statistics: true                 # Include statistical analysis
  include_convergence: true                # Include convergence analysis
  include_per_residue: true                # Include per-residue analysis
  include_hot_spots: true                  # Include hot spot analysis
  report_format: "txt"                     # txt, md, html, pdf
  include_timestamp: true                  # Include timestamp
  include_config_summary: true             # Include configuration summary
```

### 9. Debug Settings (`debug_settings`)
```yaml
debug_settings:
  debug_mode: false                        # Enable debug mode
  save_intermediate_systems: false         # Save intermediate systems
  save_energy_components: false            # Save energy components
  save_forcefield_params: false            # Save forcefield parameters
  verbose_errors: false                    # Verbose error reporting
  save_snapshots: false                    # Save trajectory snapshots
```

### 10. Reproducibility Settings (`reproducibility_settings`)
```yaml
reproducibility_settings:
  global_random_seed: 42                   # Global random seed
  save_configuration: true                 # Save exact configuration
  save_versions: true                      # Save version information
  save_system_hashes: true                 # Save system hashes
  strict_reproducibility: true             # Enable strict reproducibility
  save_environment: true                   # Save environment information
```

### 11. Performance Settings (`performance_settings`)
```yaml
performance_settings:
  memory_limit: null                       # Memory limit (GB)
  cpu_affinity: null                       # CPU affinity (list of cores)
  gpu_memory_fraction: 0.8                 # GPU memory fraction
  optimize_memory: true                    # Enable memory optimization
  batch_size: 10                           # Batch size for processing
  track_progress: true                     # Enable progress tracking
```

## 📊 Output Structure

```
mmgbsa_results/
└── complete_analysis_20250719_234006/
    ├── plots/                    # Generated plots
    ├── data/                     # Results data files
    ├── logs/                     # Log files and environment info
    ├── reports/                  # Analysis reports
    ├── cache/                    # Cache files
    └── debug/                    # Debug files (if enabled)
```

## 🔧 Common Configurations

### Quick Test
```yaml
analysis_settings:
  max_frames: 10
  run_entropy_analysis: false
  run_per_residue_decomposition: false
  use_cache: true
  verbose: 2
```

### Production Run
```yaml
analysis_settings:
  max_frames: 100
  frame_start: 1000
  frame_stride: 5
  run_entropy_analysis: true
  run_per_residue_decomposition: true
  decomp_frames: 20
  parallel_processing: true
  use_cache: true
```

### High-Precision Analysis
```yaml
analysis_settings:
  max_frames: 200
  frame_stride: 2
  run_entropy_analysis: true
  run_per_residue_decomposition: true
  decomp_frames: 50
  energy_decomposition: true

advanced_settings:
  minimization_tolerance: 1e-7
  bootstrap_samples: 5000
  nma_quality_threshold: "Excellent"
```

## 🧪 Testing

```bash
# Test configuration validation
python test_complete_config.py

# Test output structure
python complete_mmgbsa_runner.py test_complete_config.yaml

# Test with different parameters
python complete_mmgbsa_runner.py --create-config --config-name test_config.yaml
```

## 🔍 Troubleshooting

### Configuration Issues
1. **Missing sections**: Use `--create-config` to generate complete template
2. **Invalid parameters**: Check alternatives in comments
3. **File not found**: Verify input file paths

### Performance Issues
1. **Memory problems**: Reduce `max_frames` or `batch_size`
2. **Slow execution**: Enable `parallel_processing` and `use_cache`
3. **GPU issues**: Set `preferred_platform: "CPU"`

### Reproducibility Issues
1. **Different results**: Check `global_random_seed` and `deterministic_forces`
2. **Missing files**: Enable `save_configuration` and `save_environment`
3. **Version conflicts**: Check `save_versions` output

## 📚 Best Practices

1. **Always use complete configuration** for production runs
2. **Save configuration copies** for reproducibility
3. **Test with small datasets** before large runs
4. **Monitor memory usage** for large systems
5. **Use caching** for repeated analyses
6. **Validate results** before final conclusions
7. **Document parameter choices** in reports

## 🎯 Key Features

- ✅ **Complete parameter control** with alternatives
- ✅ **Reproducible analysis** with environment tracking
- ✅ **Flexible output options** in multiple formats
- ✅ **Comprehensive validation** at all stages
- ✅ **Performance optimization** settings
- ✅ **Debug and testing** capabilities
- ✅ **Professional reporting** system
- ✅ **Easy configuration** management 
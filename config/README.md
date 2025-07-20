# Configuration Files

This directory contains all configuration files for the MM/GBSA Analysis Package.

## 📁 Files

### `mmgbsa_config.yaml`
**Basic configuration file** for simple MM/GBSA analysis.
- Minimal settings for quick analysis
- Good for beginners
- Essential parameters only

### `complete_mmgbsa_config.yaml`
**Complete configuration file** with all available options.
- Full control over all parameters
- Advanced features enabled
- Detailed comments and alternatives
- Recommended for production use

## 🚀 Quick Start

### 1. Basic Analysis
```bash
# Use basic configuration
python mmgbsa_runner.py config/mmgbsa_config.yaml
```

### 2. Complete Analysis
```bash
# Use complete configuration
python complete_mmgbsa_runner.py config/complete_mmgbsa_config.yaml
```

### 3. Create New Configuration
```bash
# Create basic config
python mmgbsa_runner.py --create-config

# Create complete config
python complete_mmgbsa_runner.py --create-config
```

## 📋 Configuration Types

### Basic Configuration (`mmgbsa_config.yaml`)
```yaml
analysis_settings:
  temperature: 300
  gb_model: "OBC2"
  salt_concentration: 0.15
  max_frames: 50
  use_cache: true
```

### Complete Configuration (`complete_mmgbsa_config.yaml`)
```yaml
analysis_settings:
  # Core parameters
  temperature: 300
  gb_model: "OBC2"
  salt_concentration: 0.15
  
  # Frame selection
  max_frames: 50
  frame_start: null
  frame_end: null
  frame_stride: null
  frame_selection: "sequential"
  
  # Advanced features
  run_entropy_analysis: true
  run_per_residue_decomposition: true
  decomp_frames: 10
  energy_decomposition: false
  parallel_processing: true
```

## 🔧 Customization

### 1. Copy and Modify
```bash
# Copy basic config
cp config/mmgbsa_config.yaml my_analysis.yaml

# Edit for your needs
nano my_analysis.yaml
```

### 2. Use as Template
```bash
# Use complete config as template
cp config/complete_mmgbsa_config.yaml production_analysis.yaml

# Modify parameters
nano production_analysis.yaml
```

## 📊 Configuration Examples

### Quick Analysis
```yaml
analysis_settings:
  max_frames: 20
  frame_start: 500
  frame_stride: 10
  run_entropy_analysis: false
  run_per_residue_decomposition: false
```

### Production Analysis
```yaml
analysis_settings:
  max_frames: 100
  frame_start: 1000
  frame_stride: 5
  run_entropy_analysis: true
  run_per_residue_decomposition: true
  decomp_frames: 20
  parallel_processing: true
```

### Research Analysis
```yaml
analysis_settings:
  max_frames: 200
  frame_selection: "equidistant"
  run_entropy_analysis: true
  run_per_residue_decomposition: true
  decomp_frames: 50
  energy_decomposition: true
```

## ⚙️ Parameter Reference

### Core Parameters
- `temperature`: Analysis temperature (K)
- `gb_model`: GB model (OBC2, OBC1, HCT, GBn, GBn2)
- `salt_concentration`: Salt concentration (M)

### Frame Selection
- `max_frames`: Maximum frames to analyze
- `frame_start`: Start frame (0-indexed)
- `frame_end`: End frame (0-indexed)
- `frame_stride`: Every Nth frame
- `frame_selection`: Selection method

### Advanced Features
- `run_entropy_analysis`: Enable entropy calculation
- `run_per_residue_decomposition`: Enable residue analysis
- `decomp_frames`: Frames for decomposition
- `energy_decomposition`: Energy component analysis
- `parallel_processing`: Enable parallel processing

## 📝 Notes

- Configuration files use YAML format
- All parameters have default values
- Comments explain each parameter
- Complete config includes alternatives
- Files are validated before use

## 🔍 Validation

Configuration files are automatically validated:
- Required parameters checked
- Parameter ranges verified
- File paths validated
- Dependencies confirmed

## 📚 Documentation

For detailed configuration information, see:
- `COMPLETE_CONFIG_GUIDE.md` - Comprehensive guide
- `README.md` - Main documentation
- `CONTRIBUTING.md` - Development guide 
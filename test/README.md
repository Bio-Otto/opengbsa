# Test Directory

This directory contains all test files and configurations.

## 📁 File Structure

### 🧪 Test Scripts

- **`test_complete_config.py`**: Complete configuration test
- **`test_warning_suppression.py`**: Warning suppression test
- **`test_runner_simple.py`**: Simple runner test
- **`test_frame_selection.py`**: Frame selection test
- **`test_yaml_config.py`**: YAML configuration test

### ⚙️ Test Configurations

- **`test_complete_config.yaml`**: Complete test configuration
- **`frame_test_config.yaml`**: Frame selection test configuration
- **`frame_test_configs.yaml`**: Multiple frame selection strategies
- **`advanced_test_config.yaml`**: Advanced test configuration
- **`test_simple_config.yaml`**: Simple test configuration

### 📊 Test Data

- **`complex.pdb`**: Protein-ligand complex structure
- **`ligand.pdb`**: Isolated ligand structure
- **`ligand.sdf`**: Ligand SDF file
- **`complex.xtc`**: XTC format trajectory
- **`complex.dcd`**: DCD format trajectory

## 🚀 Running Tests

### From Main Directory

```bash
# Complete configuration test
python test/test_complete_config.py

# Warning suppression test
python test/test_warning_suppression.py

# Frame selection test
python test/test_frame_selection.py

# Simple runner test
python test/test_runner_simple.py test/frame_test_config.yaml

# YAML configuration test
python test/test_yaml_config.py
```

### From Test Directory

```bash
cd test

# Complete configuration test
python test_complete_config.py

# Warning suppression test
python test_warning_suppression.py

# Frame selection test
python test_frame_selection.py
```

## 📋 Test Configurations

### 1. Complete Test (`test_complete_config.yaml`)
```bash
python complete_mmgbsa_runner.py test/test_complete_config.yaml
```

### 2. Frame Selection Test (`frame_test_config.yaml`)
```bash
python mmgbsa_runner.py test/frame_test_config.yaml
```

### 3. Multiple Frame Strategies (`frame_test_configs.yaml`)
```bash
python test/test_frame_selection.py
```

### 4. Advanced Test (`advanced_test_config.yaml`)
```bash
python mmgbsa_runner.py test/advanced_test_config.yaml
```

### 5. Simple Test (`test_simple_config.yaml`)
```bash
python test_runner_simple.py test/test_simple_config.yaml
```

## 🧪 Test Scenarios

### 1. Configuration Validation
```bash
python test/test_complete_config.py
```
- Tests all configuration sections
- Validates parameters
- Checks output structure

### 2. Frame Selection Tests
```bash
python test/test_frame_selection.py
```
- Sequential frame selection
- Equidistant frame selection
- Random frame selection
- Frame range and stride tests

### 3. Warning Suppression
```bash
python test/test_warning_suppression.py
```
- OpenEye Toolkit warnings
- OpenMM deprecation warnings
- Import fallback tests

### 4. Simple Runner Test
```bash
python test/test_runner_simple.py test/frame_test_config.yaml
```
- YAML configuration loading
- Input file validation
- Output directory creation
- Mock analysis execution

### 5. YAML Configuration Test
```bash
python test/test_yaml_config.py
```
- YAML file loading
- Configuration validation
- Parameter checking

## 📊 Test Results

All tests should pass:
- ✅ Configuration loading
- ✅ Parameter validation
- ✅ Output structure creation
- ✅ Warning suppression
- ✅ Frame selection logic

## 🔧 Test Development

To add new tests:

1. Add test file to `test/` directory
2. Add test configuration to `test/` directory
3. Update this README
4. Make test runnable from main directory

## 📝 Notes

- All test files are organized in the `test/` directory
- Test configurations use relative paths
- Test data includes real files (complex.pdb, ligand.sdf, etc.)
- Trajectory files are large, use carefully 
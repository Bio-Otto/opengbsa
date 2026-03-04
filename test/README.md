# Test Directory

This directory contains validation scripts, unit tests, diagnostics tools, and test data for the MM/GBSA package.

## 📁 File Structure

### `test/unit/`
Unit tests for individual components.
- `test_runner_simple.py`: Basic runner tests
- `test_complete_config.py`: Configuration validation
- `test_frame_selection.py`: Frame slicing logic check
- `test_gasteiger.py`: Charge assignment validation
- `test_protein_param.py`: Protein parameterization checks

### `test/manual/`
Manual verification scripts and one-off report/visualization checks.
- `run_viz_only.py`: Generate visualization artifacts from existing results
- `check_entropy_values.py`: Quick entropy sanity checks from CSV
- `verify_interactive_report_fix.py`: Validate interactive report behavior
- `test_platform_settings.py`: Platform-setting behavior smoke script

### `test/diagnostics/`
Debugging and system capability check scripts.
- `check_charges.py`: Verify charge assignment
- `check_pbc_integrity.py`: Check simulation box
- `diagnose_energy.py`: Deep dive into energy components
- `verify_sasa.py`: Verify surface area calculations

### `test/utils/`
Helper scripts for fixing and preparing files.
- `convert_ligand.py`: Convert between formats
- `fix_pdb.py`: Fix common PDB issues
- `strip_trajectory.py`: Remove solvent/ions
- `utils/repair_trajectory.py`: Unwrap/Image trajectory

### `test/configs/`
YAML configuration files for various test cases.
- `test_complete_config.yaml`
- `frame_test_config.yaml`
- `comprehensive/*.yaml`: Comprehensive regression suite configs

### `test/data/`
Raw data files (PDB, MOL2, etc.) used by diagnostics.
(Note: Full test cases like `7khz_monomer_test` and `6xj3_test` remain as directories in `test/`)

### `test/run_analysis_test.py`
Primary manual analysis runner script.
**Usage:**
```bash
python test/run_analysis_test.py test/configs/test_complete_config.yaml
```

## 🚀 Running Tests

### Running Unit Tests
We recommend using `pytest` (if available) or running scripts individually.

```bash
# Example: Run frame selection test
python test/unit/test_frame_selection.py
```

### Running Diagnostics
Use these if you encounter specific errors.

```bash
# Check if your system has correct charges
python test/diagnostics/check_charges.py
```

### Running Comprehensive Regression
```bash
python test/run_comprehensive_tests.py
```

## 📝 Notes
- When running scripts from `test/` subdirectories, ensure the project root is in your `PYTHONPATH` or run them from the project root directory.

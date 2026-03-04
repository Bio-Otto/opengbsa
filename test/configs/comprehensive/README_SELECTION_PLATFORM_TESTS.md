# Comprehensive Tests: Selection & Platform Configuration

This document describes the new comprehensive tests added for receptor/ligand selection and platform-specific configurations.

## Test Categories

### 1. Receptor Selection Tests

#### `test_receptor_selection_chainid.yaml`
- **Purpose:** Test filtering receptor to a specific chain
- **Configuration:** `receptor_selection: "chainid 0"` (chain A only)
- **Expected Outcome:** Calculations use only the specified chain
- **Use Case:** Multi-chain proteins where you want to isolate specific chains

### 2. Ligand Selection Tests

#### `test_ligand_selection_residue.yaml`
- **Purpose:** Test custom ligand atom selection by residue
- **Configuration:** `ligand_selection: "residue 123"`
- **Additional:** Uses `"protein and not type H"` for receptor to exclude hydrogens
- **Expected Outcome:** Ligand energies computed only for specified residue atoms
- **Use Case:** Complex with multiple ligands or ligand fragments

### 3. Platform Configuration Tests

#### `test_platform_cpu.yaml`
- **Purpose:** Verify CPU-based calculations work correctly
- **Configuration:** `preferred_platform: CPU`
- **Expected Outcome:** Main analysis runs on CPU (fallback platform)
- **Use Case:** Systems without GPU, testing CPU compatibility

#### `test_platform_cpu_parallel_decomp.yaml`
- **Purpose:** Test CPU platform with parallel decomposition enabled
- **Configuration:**
  - `preferred_platform: CPU`
  - `decomposition_platform: CPU`
- **Expected Outcome:** Decomposition uses multiprocessing (safe with CPU)
- **Key Feature:** Shows logging message "✓ Running parallel decomposition on N cores"
- **Use Case:** Fast CPU-based decomposition using all available cores

#### `test_platform_cuda_cpu_decomp.yaml`
- **Purpose:** Test hybrid platform separation (GPU analysis + CPU decomposition)
- **Configuration:**
  - `preferred_platform: CUDA`
  - `decomposition_platform: CPU`
- **Expected Outcome:**
  - Main MM/GBSA analysis uses CUDA GPU acceleration
  - Decomposition uses CPU with safe multiprocessing
- **Key Feature:** User's requested feature - best performance with GPU-accelerated analysis while enabling parallelized decomposition
- **Use Case:** High-performance configurations with GPU support

### 4. Combined Tests

#### `test_selections_platform_cpu_parallel.yaml`
- **Purpose:** Integration test combining custom selections with platform parallelization
- **Configuration:**
  - `receptor_selection: "protein and not type H"`
  - `ligand_selection: "resname LIG"`
  - `preferred_platform: CPU`
  - `decomposition_platform: CPU`
  - `report_raw_energies: true`
- **Expected Outcome:**
  - Custom selections applied correctly
  - Parallel decomposition enabled
  - Raw energy values (complex, receptor, ligand) reported
- **Use Case:** Real-world scenarios requiring both custom selection and performance tuning

## How to Run Tests

### Run all comprehensive tests:
```bash
python test/run_comprehensive_tests.py
```

### Run a specific test by name:
```bash
python test/run_comprehensive_tests.py test_platform_cpu
```

### Run tests matching a pattern:
```bash
python test/run_comprehensive_tests.py test_platform_
```

## Expected Logging Output

### CPU Parallel Decomposition Enabled:
```
✓ Running parallel decomposition on 4 cores (CPU platform)
```

### Serial Decomposition (GPU Platform):
```
• Serial decomposition (decomposition_platform=CUDA incompatible with GPU-unsafe parallelization)
```

### Serial Decomposition (No Platform Set):
```
• Serial decomposition (no explicit CPU platform set; GPU parallelization disabled for safety)
```

## Platform Separation Feature Benefits

1. **Performance:** GPU acceleration for main energy calculations
2. **Parallelization:** Safe CPU-based multiprocessing for decomposition
3. **Safety:** GPU contexts isolated per process, no sharing issues
4. **Flexibility:** Users can match platform to their hardware

## Test Results Location

All test results are saved to `test/results/comprehensive/[test_name]/` with:
- Analysis results (energies, decomposition)
- Configuration copy (`analysis_config.yaml`)
- Environment info
- Logs with platform selection details

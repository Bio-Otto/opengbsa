# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Fixed
- **All 5 GB models (`HCT`/`OBC1`/`OBC2`/`GBn`/`GBn2`) are now independently validated against Amber MMPBSA.py**, not just `OBC2` (the only model prior validation work exercised). This surfaced four real, previously-unknown bugs in `mmgbsa_core.py`'s GB-model handling:
  - `GBSACalculator._create_fallback_obc_force` always built an `OBC2`-physics force regardless of the configured `gb_model` -- `HCT`/`OBC1`/`GBn`/`GBn2` silently fell back to `OBC2` energetics whenever this fallback path was taken (which, for Coordinate Mode, was every run, since `openmmforcefields`'s `SystemGenerator` always rejects the `implicitSolvent` kwarg there). Fixed by dispatching to the matching `openmm.app.internal.customgbforces` class per model.
  - Explicit `receptor_topology`/`ligand_topology` inputs (loaded from separate prmtop files in Native Mode) never had their GB radii synchronized with the complex structure's radii, producing a complex/receptor+ligand radii mismatch large enough to change a real system's `OBC2` result by ~200 kcal/mol.
  - `GBn`'s radii set was incorrectly mapped to `mbondi` instead of `bondi` -- confirmed this wasn't OpenGBSA-specific: ParmEd's own official `createSystem(implicitSolvent=app.GBn)` path raises the same "Radii must be between 1 and 2 Angstroms for neck lookup" error when given mbondi's valid Amber-standard hydroxyl-hydrogen radii, since `GBn`'s neck-lookup table is calibrated for `bondi`'s radius range specifically.
  - `GBn`/`GBn2` need a model-specific per-element `screen` (Born-radius scaling) value from OpenMM's own lookup table, not the generic prmtop-derived value every other model uses -- this was the largest remaining error source, responsible for a ~13 kcal/mol `GBn2` discrepancy against Amber before being fixed.
  See `test/configs/ppi_1gcq_test/README.md`'s "All 5 GB models, validated" section for the full comparison table and per-bug detail.

### Added
- `test/configs/ppi_1gcq_test/prepare_system.py`: new `--gb-model {HCT,OBC1,OBC2,GBn,GBn2}` flag for `--build-amber-reference`, so the independent Amber reference topology's radii set can be built matching any of the 5 supported GB models (previously hardcoded to whatever `tleap`/ParmEd defaulted to, implicitly OBC2-only).

## [0.0.7] - 2026-09-03

### Added
- **Nucleic-acid support**: `binding_mode: ppi` and Native Mode topology-splitting now correctly handle RNA and DNA "ligands" (previously, ligand-detection logic assumed a small-molecule `ligand_resname` and silently mis-split protein-RNA/DNA complexes). Validated end-to-end against independent Amber references on real protein-RNA and protein-DNA datasets.
- Five new self-contained, independently-verified validation test cases under `test/configs/` (NAMD/CHARMM peptide, GROMACS PPI x2, RNA-protein, DNA-protein), each with its own data-fetch script and expected-result README, plus a new ReadTheDocs Tutorial section walking through all five.
- `opengbsa --fetch-test-data`/`--list-test-data` CLI flags (Pooch + Zenodo-backed) so `test/unit`'s core fixtures no longer need to be committed to the repository. Published as Zenodo record [10.5281/zenodo.22139199](https://doi.org/10.5281/zenodo.22139199).
- `CITATION.cff` for software citation metadata.

### Changed
- **Performance**: the main per-frame MM/GBSA energy loop is now parallelized across CPU cores (`ProcessPoolExecutor`), and redundant duplicate energy queries were removed from the per-frame calculation (~42 -> ~24 `getState()` calls/frame). Measured ~2.37x speedup on a 300-frame run, with numerically-verified identical results (bit-identical VdW/electrostatic/surface-area/internal energies; GB shows ~1e-4 kcal/mol differences that are OpenMM's own floating-point non-determinism, present even between two runs of the unmodified original code). Only activates on the CPU platform when entropy_method is not `quasiharmonic`; GPU platforms and QHA runs are unaffected.
- Documented the previously-undocumented `nonbonded_cutoff` config option in the configuration reference (no default behavior change).
- `test/unit/`: added 39 new tests covering `ConfigManager`, `TopologyLoader`'s Amber/GROMACS loading paths, and an end-to-end `MMGBSARunner` regression test (18 -> 57 tests). Moved three files that contained no actual test functions (`test_protein_param.py`, `test_numpy_patch.py`, `test_runner_simple.py` -- developer debug scripts pytest was silently collecting without running anything) to `test/manual/`.

### Fixed
- **Removed dead `mmgbsa/core.py`:** the old monolithic module was fully shadowed by the `mmgbsa/core/` package and unreachable; deleted (~4,570 lines).
- **Removed dead `mmgbsa/core/` package** (not the same file as the entry above -- this is the *package* that had shadowed it): despite its docstrings describing a "modular refactor" of `mmgbsa_core.py` into `platform`/`caching`/`topology`/`parameterization`/`analysis`/`results`/`calculator` submodules, the package's own code admitted this refactor was never actually wired in (`AnalysisEngine`'s `frame_selector`/`energy_calc`/`decomposition`/`aggregator` were built but never called; its `GBSACalculator` was a thin facade whose every method delegated straight back to the real `mmgbsa.mmgbsa_core.GBSACalculator`). No production code outside `mmgbsa/core/` itself used any of its manager classes -- only a handful of `scripts/scratch/` developer scripts (also removed) did. `mmgbsa/runner.py`, `mmgbsa/complete_runner.py`, and `mmgbsa/decomposition.py` now import `GBSACalculator` directly from `mmgbsa.mmgbsa_core`, with no API or behavior change (the facade's `__init__` signature was already an exact passthrough copy of the real one).
- **Native TPR mode (`tpr_loader.py`):** bonded parameters (bonds/angles/dihedrals/impropers) are now parsed from the `.tpr` file instead of being fabricated; triclinic box vectors are now computed correctly. CMAP terms are detected and warned about (not yet supported for native TPR).
- **GB-derived per-residue decomposition:** implemented a full analytical OBC2 decomposition (`_gb_solvation_decomposition`) and fixed a ~50% systematic undercounting bug in the Gohlke-Kollman pair-energy split when the decomposition target (e.g. the ligand) has no entry of its own in the residue map.
- **vdW decomposition:** the extracted CHARMM NBFIX `acoef`/`bcoef` pair table is now actually used in the standalone per-residue pairwise vdW calculation (previously silently fell back to Lorentz-Berthelot combining rules even when a real NBFIX table was available).
- **Electrostatic decomposition:** removed an incorrectly-applied Debye-Huckel salt-screening factor from the pairwise electrostatic decomposition term.
- **Salt screening (kappa):** `OBC1`/`OBC2` GB models now correctly apply Debye-Huckel salt screening when `salt_concentration > 0`, via a `CustomGBForce`-based `GBSAOBC2Force`; previously this was silently ignored for every run.
- **GB radii for GROMACS-origin systems:** the configured GB radii set (e.g. `mbondi2` for OBC2) is now actually applied when converting a GROMACS `.top` to Amber format (`GromacsPreprocessor.convert_to_amber`), which is the code path real GROMACS-input runs go through.
- **1-4 nonbonded exception detection:** fixed a heuristic that only sampled the first 5 `NonbondedForce` exceptions to decide whether Amber-style 1-4 scaling needed to be substituted in; it now scans all exceptions, and prefers the structure's own CHARMM-native `adjusts` list over hardcoded Amber SCEE/SCNB values when available.
- **Force-group assignment:** fixed a collision where all `CustomBondForce` instances (including the 1-4 VDW correction) were unconditionally reassigned to one force group, and where CHARMM improper torsions collided with the surface-area force's group; both are now identified by their energy-function signature and given dedicated groups.
- **Interaction entropy:** rewrote using `scipy.special.logsumexp` for numerical stability; now raises instead of silently returning 0.0 on failure, and warns when `sigma(dE)` exceeds the reliability threshold from Duan et al. 2016.
- **Normal-mode vibrational entropy (`getVibrationalEntropyCM`):** corrected a formula that omitted `hbar` entirely.
- **Output directory resolution:** fixed a `KeyError` crash in `runner.py`/`complete_runner.py` when a config used the documented `params.output_directory` key instead of an undocumented alternate location.
- Removed two non-importable orphaned code fragments (`mmgbsa/combined_system.py`, `mmgbsa/config_patch.py`) left over from the `mmgbsa/core/` package extraction; their logic already exists properly in `mmgbsa_core.py`/`config.py`.
- Removed stale root-level scratch scripts (`run_mmpbsa.py`, `mmgbsa_cli.py`) that imported the deleted `mmgbsa.core` module and duplicated the packaged `opengbsa` CLI entry point; installation docs now reference `opengbsa --help` directly.
- Removed two broken `test/unit/` scripts (`test_combined_system.py`, `test_gasteiger.py`) that referenced a hardcoded developer machine path and the deleted `mmgbsa.core` module, breaking `pytest test/` collection for anyone running it from a fresh checkout.

### Known issues
- `TopologyLoader._load_gromacs`'s `.tpr` branch (`mmgbsa/topology.py`) calls `parmed.load_file()` directly on a `.tpr` file, but ParmEd has no `.tpr` parser and raises `FormatNotFound`. This code path is separate from (and apparently unreachable compared to) the real native-TPR support `mmgbsa_core.py` actually uses (`tpr_loader.py`'s `TprParser`-based loader). Documented in `test/unit/test_topology_loading.py`; not yet fixed or removed.

## [0.0.6] - 2026-03-04

### Added
- **Comprehensive Test Suite:** 48 automated YAML-based test configurations covering all major features
  - All GB models (OBC1, OBC2, HCT, GBn, GBn2), SA models (ACE, LCPO), entropy methods, decomposition variants
  - Platform-specific tests: CPU, CPU parallel decomposition, CUDA+CPU, OpenCL+CPU
  - Ligand/receptor selection tests, dimer mode, DCD/GRO/PRMTOP/TPR input format tests
  - Automated test runner (`test/run_comprehensive_tests.py`) with early-exit on ≥3 failures
- **Per-Residue Decomposition Platform Settings:** All test configs now use `decomposition_platform: CPU` with `parallel_processing: true` for multiprocess decomposition
- **GPU Platform Randomization:** Comprehensive configs randomize `preferred_platform` between CUDA and OpenCL
- **Dimer Mode Test:** `dimer_dual_test.yaml` validates protein-protein-ligand dimer analysis end-to-end
- **`scripts/update_configs.py`:** Utility script to batch-update test config settings (decomp_frames, platforms)

### Fixed
- **All 48 test configs:** Corrected placeholder paths (`test/complex.pdb`, `test/complex.dcd`) to real test data
- **Integer temperature validation:** Changed all integer `temperature: 300` to `temperature: 300.0` (float required)
- **`decomp_frames > max_frames` errors:** Reduced decomp_frames to be ≤ max_frames across advanced configs
- **`dimer_dual_test`:** Fixed "Ligand residues not found" by adding `ligand_resname: LIG` and switching from solvated PDB (118K lines) to `.tpr` input for automatic solvent stripping
- **`frame_test_configs`:** Rewrote non-standard multi-section YAML as a flat standard config
- **`6xj3_config` family:** Corrected paths from `analysis_6xj3/` to `test/data/6xj3_pdb_test/`
- **`advanced_test_config`:** Fixed file paths and integer temperature

### Changed
- **Test configs consolidated:** All test YAML files moved from `test/configs/` to `test/configs/comprehensive/`
- **README.md:** Updated to reflect current CLI (`mmgbsa` / `python -m mmgbsa.cli`), platform settings, test suite, and expanded error table
- **`test/run_comprehensive_tests.py`:** Now requires the `opengbsa` conda environment for openff-toolkit support

## [0.0.5] - 2026-01-26


### Added
- **LCPO Surface Area Model**: Implemented Linear Combinations of Pairwise Overlaps method
  - More accurate physics-based surface area calculations
  - Compatible with AMBER MMPBSA.py default behavior
  - Configurable via `sa_model: LCPO` in analysis settings
- **OpenMM Source Build Integration**: Added support for OpenMM 8.1+ with LCPO
- **Automatic Parameter Assignment**: Uses `getLCPOParamsTopology()` for robust atom typing
- **Force Group Energy Extraction**: Proper SA energy isolation for LCPO calculations

### Fixed
- **LCPO Zero Energy Bug**: Critical fixes for LCPO implementation
  - Added probe radius (1.4 Å) to particle radii (matching reference implementation)
  - Added force group assignment (`setForceGroup(4)`) for proper energy extraction
  - Verified against standalone tests (45.62 kcal/mol) and integrated GBSA (65.69 kcal/mol)
- **Frame-by-frame decomposition CSV**: Fixed missing CSV file generation for HTML reports
- **VdW energy extraction**: Fixed parameter injection for decomposition module

### Changed
- **Surface Area Model Selection**: Users can now choose between ACE (fast) and LCPO (accurate)
- **Documentation**: Updated README, CONFIGURATION, and added LCPO usage guides

### Performance
- **LCPO vs ACE Comparison** (6T1H system, Frame 0):
  - ACE: SA Complex=124.40 kcal/mol, Delta=-13.92 kcal/mol
  - LCPO: SA Complex=65.69 kcal/mol, Delta=-8.17 kcal/mol
  - Computation time similar for single-frame analysis


### Added
- **Complete MM/GBSA Analysis Package** with advanced features
- **Multiple GB Models**: OBC2, OBC1, HCT, GBn, GBn2 support
- **Normal Mode Analysis**: Entropy calculations with ultra-robust minimization
- **Per-Residue Decomposition**: Detailed residue-ligand interaction analysis
- **YAML Configuration**: Single configuration file for all parameters
- **Advanced Validation**: Input validation and result quality checks
- **Parallel Processing**: Multi-core support for faster analysis
- **Caching System**: Reuse prepared systems for efficiency
- **Comprehensive Reporting**: Detailed analysis reports with plots
- **Docker Support**: Containerized deployment
- **Frame Selection Strategies**: Sequential, equidistant, and random sampling
- **Energy Decomposition**: Component-wise energy analysis
- **Hot Spot Identification**: Key binding site analysis
- **Advanced Visualization**: Publication-quality plots and charts
- **ProLIF Integration**: Protein-ligand interaction fingerprinting
- **GitHub Issue #185 Solution**: HTML-based interaction network export

### Changed
- **Complete rewrite** of MM/GBSA analysis pipeline
- **Enhanced error handling** and validation
- **Improved performance** with parallel processing
- **Better documentation** with comprehensive guides

### Fixed
- **Memory optimization** for large trajectories
- **Numerical stability** in entropy calculations
- **File handling** issues with various formats
- **Configuration validation** and error reporting

### Technical Details
- **Python 3.8+** compatibility
- **OpenMM 8.0+** integration
- **CUDA support** for GPU acceleration
- **Cross-platform** compatibility
- **MIT License** for open source use

## [0.0.3] - 2024-01-15

### Added
- Initial MM/GBSA implementation
- Basic trajectory analysis
- Simple energy calculations

### Changed
- Improved code structure
- Better error handling

## [0.0.2] - 2024-01-10

### Added
- Basic molecular dynamics analysis
- Trajectory processing capabilities

## [0.0.1] - 2024-01-05

### Added
- Initial project setup
- Basic file structure
- README documentation

---

## Version History

- **v0.0.4**: Complete MM/GBSA analysis package with advanced features
- **v0.0.3**: Initial MM/GBSA implementation
- **v0.0.2**: Basic molecular dynamics analysis
- **v0.0.1**: Project initialization

## Future Plans

### v0.1.0 (Planned)
- **Enhanced GUI**: Web-based interface
- **More Force Fields**: Additional force field support
- **Cloud Integration**: AWS/Azure deployment options
- **API Development**: RESTful API for programmatic access

### v0.2.0 (Planned)
- **Machine Learning**: ML-based binding affinity prediction
- **Advanced Analytics**: Statistical analysis tools
- **Plugin System**: Extensible architecture
- **Performance Optimization**: Further speed improvements 
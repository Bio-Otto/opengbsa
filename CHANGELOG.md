# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.0.4] - 2024-01-20
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
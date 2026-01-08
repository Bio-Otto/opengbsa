# Changelog

All notable changes to this project will be documented in this file.

## [v0.0.5] - 2026-01-08

### 🚀 Major Features
*   **Multi-Engine Support:**
    *   **OpenMM Native:** Added support for loading serialized OpenMM XML Systems.
    *   **Gromacs Integration:** Full automatic conversion pipeline (.top -> .prmtop) integrated.
    *   **Amber Native:** Improved InputManager to robustly handle Amber prmtop/inpcrd files.
    *   **Documentation:** Created comprehensive MULTI_ENGINE_GUIDE.md.

### 📊 Visualization Enhancements
*   **Statistical Plots:**
    *   rolling_average: Monitors energy stabilization over time.
    *   energy_convergence: Cumulative average analysis.
    *   components_pie_chart: Visual breakdown of dG components.
    *   entropy_convergence: Interaction Entropy stability tracking.
*   **HTML Reporting:**
    *   Integrated all new plots into the advanced_analysis_report.html.

### 📚 Documentation (Internationalization)
*   **English Translation:** All documentation and comments are fully in English.
*   **New Guides:** INSTALLATION.md, CONFIGURATION.md, GROMACS_GUIDE.md, OUTPUTS.md.
*   **Read The Docs:** Updated conf.py to support Markdown and fixed Math rendering.

### 🐛 Bug Fixes
*   Fixed zero VdW values in Per-Residue Decomposition.
*   Fixed atom count mismatch during visualization frame extraction.

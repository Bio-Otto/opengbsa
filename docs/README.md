# OpenGBSA Documentation

Welcome to the comprehensive user guide for OpenGBSA.

This folder contains specialized guides for different modules and usage scenarios of the project.

---

## 📚 Guides

### 🚀 Getting Started
*   **[Installation Guide](INSTALLATION.md):**
    *   Setting up the Conda environment, installing required libraries (OpenMM, MDTraj), and verification steps.

### ⚙️ Settings and Configuration
*   **[Configuration Guide](CONFIGURATION.md):**
    *   Structure of `config.yaml`, detailed explanations of `input`, `analysis_settings`, and `visualization` parameters.
*   **[GROMACS User Guide](GROMACS_GUIDE.md):**
    *   Special considerations when working with GROMACS `.top` and `.xtc` files (Handling .itp files, atom counts, etc.).
*   **[Multi-Engine Support (AMBER, NAMD, OpenMM)](MULTI_ENGINE_GUIDE.md):**
    *   Instructions for using other simulation engines and formats.

### 📊 Results and Analysis
*   **[Output Files and Reporting](OUTPUTS.md):**
    *   Meanings of the generated HTML report, CSV files, and plots (`rolling_average`, `convergence`, `heatmap`).
    *   How to use the PyMOL visualization file (`view_binding.pml`).

---

## 🏗️ Project Structure

*   `mmgbsa/`: Source codes.
    *   `core.py`: Main calculation engine.
    *   `visualization.py`: Charting and reporting module.
    *   `conversion.py`: Format converter (Gromacs -> Amber).
*   `run_mmpbsa.py`: Main executable script.
*   `config_master.yaml`: Master template file containing all parameters.
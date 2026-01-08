# Output Files and Reporting

When the OpenGBSA analysis is complete, several files are created in the `results` directory (or the `output_directory` specified in your config).

---

## 1. Reports

### HTML Report (`advanced_analysis_report_LIG.html`)
**The most important file.** An interactive, visually enriched summary of the analysis.
*   **Analysis Summary:** Total binding energy and statistics.
*   **Energy Statistics:** Moving average, convergence, and component distribution plots.
*   **Energy Heatmap:** A heatmap showing which amino acids contribute strongly to binding.

### CSV Results File (`fixed_enhanced_mmgbsa_results_OBC2.csv`)
Contains raw data. Can be used for detailed analysis with Excel or Pandas.
*   `frame`: Frame number analyzed.
*   `binding_energy`: Total binding energy (kcal/mol).
*   `complex_energy`, `protein_energy`, `ligand_energy`: Potential energies of components.
*   `delta_nb`, `delta_gb`, `delta_sa`: Differences (Delta) per component.

---

## 2. Plots (`plots/` Directory)

### `rolling_average_LIG.png`
**Purpose:** To understand if the system has reached equilibrium (stabilization).
*   **Gray Line:** Raw data (Frame-by-frame).
*   **Dark Line:** Moving average (Shows the trend).

### `convergence_plot_LIG.png`
**Purpose:** To prove that the simulation duration is sufficient.
*   **Green Line:** Cumulative average. If the graph flattens out (plateaus) towards the end, the result is reliable.

### `components_pie_LIG.png`
**Purpose:** To understand the source of binding energy.
*   **VDW:** Van der Waals contributions (Shape complementarity).
*   **Electrostatic:** Charge interactions.
*   **Solvation (GB/SA):** Solvent effects.

### `energy_heatmap_LIG.png`
**Purpose:** To identify critical amino acids (Hot spots).
*   **Red/Blue Colors:** Negative (red) areas represent regions binding tightly to the ligand.

---

## 3. Visualization (`view_binding.pml`)

An automatically generated session file for the PyMOL program.
When you open this file with PyMOL:
1.  Protein and Ligand are loaded automatically.
2.  Important interactions (Hydrogen bonds, etc.) are shown with yellow dashed lines.
3.  Critical amino acids are highlighted in "sticks" format.

**Usage:**
From terminal: `pymol results/view_binding.pml`

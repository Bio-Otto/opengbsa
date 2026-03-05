Visualization
=============

OpenGBSA emphasizes visual analysis to make sense of complex energy data.

3D Structure Viewer
-------------------

The interactive report includes a 3D visualization of your protein-ligand complex.

- **Color Coding**: Residues are colored by their contribution to binding.
  - **Red**: Unfavorable (repulsive) interaction.
  - **White**: Neutral.
  - **Blue/Green**: Favorable (attractive) interaction.
- **Interactivity**: Click and drag to rotate, scroll to zoom. Hover over atoms to see residue names.

Energy Plots
------------

1. **Binding Energy Time Series**: Shows stability of binding over the simulation trajectory. A stable plateau indicates convergence.
2. **Component Breakdown**: Bar chart showing the average contribution of VdW, Electrostatic, and Solvation terms. This reveals the "driving force" of binding (e.g., is it hydrophobic or electrostatic?).
3. **Residue Decomposition**: "Manhattan Plot" style bar chart. Bars extending downwards (negative) are key binding residues (Hotspots).

.. image:: ../_static/logo.png
   :align: center
   :width: 400px
   :alt: Placeholder for actual plot image

Generating Publication-Quality Figures
--------------------------------------

The static PNGs in the ``figures/`` folder are generated at high DPI (300 dpi) using Matplotlib/Seaborn styles. They are ready for direct inclusion in manuscripts or presentations.

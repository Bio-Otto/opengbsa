# -- Path setup --------------------------------------------------------------
import os
import sys
sys.path.insert(0, os.path.abspath('../..'))

# -- Project information -----------------------------------------------------
project = 'OpenGBSA'
copyright = '2025, Ozbek Lab'
author = 'H. Ibrahim Ozdemir'
release = '0.0.5'

# -- Mock Imports ------------------------------------------------------------
autodoc_mock_imports = [
    'openmm', 
    'openmm.app', 
    'openmm.unit',
    'simtk',
    'simtk.openmm',
    'simtk.unit',
    'mdtraj',
    'openff',
    'openff.toolkit',
    'rdkit',
    'rdkit.Chem',
    'prolif',
    'MDAnalysis',
    'numpy',
    'pandas',
    'scipy',
    'matplotlib',
    'seaborn'
]

# -- General configuration ---------------------------------------------------
extensions = [
    'myst_parser',
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
]

# Enable Markdown
source_suffix = {
    '.rst': 'restructuredtext',
    '.md': 'markdown',
}

# -- Options for HTML output -------------------------------------------------

# Use the Furo theme for a modern, Material-like appearance
html_theme = "furo"  # Chosen for its clean look and excellent admonition support

html_logo = "_static/logo.png"
html_static_path = ["_static"]

title = "A Fully Open-Source MM/GBSA Analysis Tool"
html_title = "A Fully Open-Source MM/GBSA Analysis Tool"
html_short_title = "OpenGBSA"
html_theme_options = {
    "light_css_variables": {},
    "dark_css_variables": {},
    "default_light_mode": True,
    "default_color_mode": "light",
}
templates_path = ["_templates"]
html_sidebars = {
    "**": [
        "sidebar/brand.html",
        "sidebar/navigation.html",
        "sidebar/variant-selector.html",
        "sidebar/search.html",
    ]
}

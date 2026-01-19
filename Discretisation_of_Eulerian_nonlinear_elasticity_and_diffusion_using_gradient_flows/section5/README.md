# Porosity waves model

This directory contains a FEniCS-based simulation for quasi-static Euler-Euler mixture model with solid and fluid momentum balance.

**Authors:** Andrea Zafferi, Dirk Peschka

**Note:** Slight code conversion and documentations were done with AI assistance using GitHub Copilot.

## Files

- **`Sim_Section5_v01.ipynb`** - Jupyter notebook containing the main simulation for the two-phase nonlinear elastic and diffusive material in Eulerian coordinates.
- **`Postprocessing_Section5_v01.ipynb`** - Jupyter notebook containing the postprocessing the two-phase nonlinear elastic and diffusive material in Eulerian coordinates.
- **`utilities.py`** - Utility functions for I/O and data handling

## Usage

### Basic Usage
Run the Jupyter notebook:

```bash
jupyter notebook Sim_Section5_v01.ipynb
```

```bash
jupyter Postprocessing_Section5_v01.ipynb
```

Run the simulation notebook *Sim_Section5_v01.ipynb* to generate a file containing the PDE parameters (saved as *.json*) and the solutions (saved as *.h5*) into *output* folder.
Run the postprocessing notebook *Postprocessing_Section5_v01.ipynb* to load the files previously generated in the *output* folder and plot time snapshot for the evolution of the fluid content, the fluid velocity and the determinant of the solid deformation gradient.


## Key Features

- **Energy formulation** for elastic-diffusive materials
- **Mixed finite element spaces** for the solid and fluid displacements, the fluid content and the chemical potential
- **Energy-consistent discretization** with proper dissipation handling

## Requirements

- FEniCS (dolfin)
- Jupyter Notebook or JupyterLab
- NumPy, Matplotlib
- tqdm (for progress bars)
- Standard Python libraries (json, os, sys, logging)
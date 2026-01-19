# Dissipative hyperelastic solid

This directory contains a FEniCS-based simulation for nonlinear elastic materials, used in.

**Authors:** Andrea Zafferi, Dirk Peschka

**Note:** Slight code conversion and documentations were done with AI assistance using GitHub Copilot.

## Files

- **`stationary_sim.ipynb`** - Jupyter notebook containing the main simulation and the postprocessing for the study of stationary solutions of the nonlinear elastic solid both in Lagrangian and Eulerian coordinates.
- **`darcy_sim.ipynb`** - Jupyter notebook containing the main simulation and the postprocessing for the study of evolutionary solutions of the nonlinear elastic solid with Darcy-like dissipation both in Lagrangian and Eulerian coordinates.
- **`utils.py`** - Utility functions for I/O and data handling


## Basic Usage
Run the Jupyter notebook:

```bash
jupyter notebook postprocessing_wave.ipynb
```

The first cell of each notebook solves the associated PDE problem. The remaining cells are for postprocessing.


## Key Features

- **Energy formulation** for elastic materials
- **Mixed finite element spaces** for displacement
- **Energy-consistent discretization** with proper dissipation handling

## Requirements

- FEniCS (dolfin)
- NumPy, Matplotlib
- tqdm (for progress bars)
- Standard Python libraries (json, os, sys, logging)
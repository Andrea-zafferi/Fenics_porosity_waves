# Porosity Wave Simulation

This directory contains a FEniCS-based simulation for elastic-diffusive material with porosity waves, converted from the original Jupyter notebook `EUL_elastic-diffusive-y.ipynb`.

**Authors:** Andrea Zafferi, Dirk Peschka

**Note:** Slight code conversion and documentations were done with AI assistance using GitHub Copilot.

## Files

- **`porosity_wave.py`** - Main simulation script (converted from notebook)
- **`postprocessing_wave.ipynb`** - Postprocessing and visualization notebook
- **`utils.py`** - Utility functions for I/O and data handling
- **`parameter_sets/`** - Directory containing parameter files
  - `pars_diffusive.json` - Parameters for diffusive regime
  - `pars_flowing.json` - Parameters for flowing regime

## Usage

### Basic Usage

Run a simulation with a specific parameter file:

```bash
python porosity_wave.py parameter_sets/pars_diffusive.json
```

or

```bash
python porosity_wave.py parameter_sets/pars_flowing.json
```

### Postprocessing and Visualization

After running simulations, use the Jupyter notebook for analysis:

```bash
jupyter notebook postprocessing_wave.ipynb
```

The postprocessing notebook provides:
- Energy evolution analysis and plotting
- Deformation and concentration field visualization  
- Dissipation rate calculations
- Comparison between different parameter sets
- Animation frame generation capabilities

## Parameter Files

The simulation uses JSON parameter files that specify:

- **Material parameters**: shear modulus (μ), diffusion coefficient (mob), viscosity (visc)
- **Physical parameters**: gravity (g0), gradient term (ε), reference concentration (c0)
- **Simulation parameters**: time steps (n_steps), final time (T), mesh resolution (k)
- **Initial conditions**: bubble parameters (α, mag)
- **Output settings**: basename and output directory (fout)

## Output

The simulation generates:

1. **Console output**: Progress information and simulation status
2. **Data files**: State files saved as HDF5 + JSON pairs in the output directory
   - Format: `out_<basename>/stateP1_<k>_<n_steps>_<step_number>.h5/.json`

## Key Features

- **Eulerian energy formulation** for elastic-diffusive materials
- **Mixed finite element spaces** for displacement, concentration, and chemical potential
- **Energy-consistent discretization** with proper dissipation handling

## Requirements

- FEniCS (dolfin)
- NumPy
- tqdm (for progress bars)
- Standard Python libraries (json, os, sys, logging)
# AITW-Darcy-Infiltration-Model

This project contains code for simulating infiltration processes in porous materials using the FEniCSx finite element framework.

## Files

- **`data.py`** - Configuration file containing input parameters for the simulation
- **`infiltration.py`** - Main simulation script implementing the infiltration model

## Description

The simulation models capillary infiltration in porous media using:
- 3D finite element mesh (hexahedral elements)
- Nonlinear saturation curve based on pore size distribution
- Adaptive timestepping for convergence
- Parallel computing support via MPI

### Key Features

- **Mesh**: 3D rectangular domain with configurable dimensions
- **Physics**: Capillary pressure, surface tension, contact angle effects
- **Material properties**: Porosity, permeability (anisotropic), pore radius distribution
- **Solver**: Newton-Raphson with adaptive timestepping
- **Output**: VTX files for visualization and numpy arrays for data analysis

## Requirements

This project requires FEniCSx and related dependencies with:
- `fenics-dolfinx` - Main FEniCS library
- `petsc4py` - Parallel linear algebra
- `mpi4py` - Message passing interface
- `matplotlib` - Plotting and visualization
- `numpy` - Numerical computing
- `tqdm` - Progress bars

## Usage

Run the simulation by providing the parameter file:

```bash
python infiltration.py data.py
```

## Output

The simulation generates:
- `{label}_solution.npy` - Complete solution array for analysis
- `{label}_output.bp` - VTX mesh files for ParaView visualization
- `{label}.log` - Simulation log file

## Parameters

Key simulation parameters (configurable in `data.py`):
- Mesh dimensions and resolution
- Material properties (porosity, permeability, surface tension)
- Fluid properties (viscosity, contact angle)
- Pore size distribution parameters
- Timestepping control
- Boundary conditions

---

*Part of the AITW project at VTT Technical Research Centre of Finland*

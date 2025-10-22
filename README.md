# AITW-Darcy-Infiltration-Model

A 3D computational model for simulating infiltration processes in porous media using Darcy's law. This model is implemented using FEniCSx (DOLFINx) and solves the nonlinear infiltration problem with adaptive timestepping.

## Overview

This project simulates fluid infiltration in porous materials by solving a coupled system of equations involving:
- Capillary pressure and saturation relationships
- Darcy flow with anisotropic permeability
- Time-dependent saturation evolution
- Contact angle effects and surface tension

The simulation uses a Newton solver with adaptive timestepping to handle the nonlinear nature of the problem and ensure numerical stability.

## Features

- **3D hexahedral mesh** for spatial discretization
- **Adaptive timestepping** for efficient convergence
- **Nonlinear saturation curves** based on capillary radius distribution
- **Anisotropic permeability** (longitudinal and transverse)
- **MPI parallelization** support for large-scale simulations
- **VTX output** for visualization in ParaView or similar tools
- **Progress monitoring** with tqdm

## Requirements

### Dependencies

- Python 3.7+
- [FEniCSx](https://fenicsproject.org/) (DOLFINx)
- petsc4py
- mpi4py
- NumPy
- matplotlib
- tqdm
- UFL (Unified Form Language)
- Basix

### Installation

Install FEniCSx following the [official installation guide](https://github.com/FEniCS/dolfinx#installation).

For other dependencies:
```bash
pip install numpy matplotlib tqdm
```

Note: `petsc4py`, `mpi4py`, `ufl`, and `basix` are typically installed as part of the FEniCSx installation.

## Usage

### Basic Usage

Run the simulation with a parameter file:

```bash
python infiltration.py data.py
```

Or with MPI for parallel execution:

```bash
mpirun -n 4 python infiltration.py data.py
```

### Parameter Configuration

The simulation parameters are defined in a Python file (e.g., `data.py`). Key parameters include:

#### Mesh Parameters
- `xmin`, `xmax`, `ymin`, `ymax`, `zmin`, `zmax`: Domain dimensions (meters)
- `nx`, `ny`, `nz`: Number of mesh elements in each direction

#### Physical Parameters
- `ϕ`: Porosity (0-1)
- `γ`: Surface tension (N/m)
- `θdeg`: Contact angle (degrees)
- `r_μ`: Mean pore radius (m)
- `rstd`: Relative standard deviation of pore radius
- `r_cutoff`: Minimum pore radius cutoff (m)
- `k_long`, `k_trans`: Longitudinal and transverse permeability (m²)
- `μ`: Fluid viscosity (Pa·s)
- `u_ext`: External pressure boundary condition (Pa)

#### Simulation Parameters
- `Δt0`: Initial timestep size (s)
- `Δt_min`, `Δt_max`: Minimum and maximum timestep bounds (s)
- `t_end`: Maximum simulation time (s)
- `label`: Simulation identifier for output files

### Example Parameter File

See `data.py` for a complete example configuration representing infiltration in a 3D plate geometry.

## Output

The simulation generates the following output files:

1. **`{label}_solution.npy`**: NumPy array containing the solution at each timestep (for serial runs only)
2. **`{label}_output.bp`**: VTX format file for visualization (can be opened in ParaView)
3. **`{label}.log`**: Detailed log file with solver information and convergence details

## Visualization

Open the `.bp` output file in ParaView to visualize:
- Pressure field (`u`)
- Saturation field (`S`)
- Time evolution of infiltration

## Mathematical Model

The model solves the following weak form:

```
∫ δu · ϕ(S - S_n) dx - Δt ∫ ∇δu · J dx = 0
```

where:
- `u`: Capillary pressure
- `S`: Saturation (function of capillary radius)
- `J`: Darcy flux: J = -(k_s · k / μ) · ∇u
- `k_s`: Saturation-dependent permeability factor
- `ϕ`: Porosity

The saturation curve is computed using an error function based on the pore radius distribution.

## Solver Settings

The Newton solver uses:
- Absolute tolerance: 1e-11
- Maximum iterations: 10
- Relaxation parameter: 0.8 (adaptive)
- Linear solver: MUMPS direct solver

## License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE](LICENSE) file for details.

## Contributing

Contributions are welcome! Please feel free to submit issues or pull requests.

## Contact

For questions or support, please open an issue in the repository.
import json
from dataclasses import dataclass


@dataclass
class Params:
    # Label for the simulation
    label: str = '3d_plate'

    # Mesh dimensions (m)
    xmin: float = 0
    xmax: float = 5e-2
    ymin: float = 0
    ymax: float = xmax
    zmin: float = 0
    zmax: float = 1e-3
    # Number of nodes
    nx: int = 40
    ny: int = 160
    nz: int = 40

    # Timestepping
    # Initial time step length (s)
    delta_t0: float = 1e-1
    delta_tmin: float = 1e-2
    delta_tmax: float = 10
    # Simulation time (s)
    t_end: float = 120

    # Porosity
    porosity: float = 0.6

    # Surface tension (N/m)
    surface_tension: float = 30e-3

    # Contact angle
    contact_angle: float = 10

    # Mean pore radius (m)
    mean_pore_radius: float = 20e-6

    # Relative standard deviation of pore radius
    rstd: float = 1/5

    # Pore radius cutoff
    r_cutoff: float = 2e-6

    # Permeability
    k_long: float = 1e-12
    k_trans: float = 1e-16

    # Fluid viscosity (Pa s)
    viscosity: float = 1e-3

    # External pressure
    u_ext: float = -1e-6

    def to_json(self, json_file: str):
        """Save the parameters to a JSON file"""
        data = {k: v for k, v in self.__dict__.items() if not k.startswith('_')}
        with open(json_file, 'w') as f:
            json.dump(data, f, indent=4)

    @classmethod
    def from_json(cls, json_file: str) -> 'Params':
        """Create an instance from a JSON file"""
        with open(json_file, 'r') as f:
            data = json.load(f)
        return cls.from_dict(data)

    @classmethod
    def from_dict(cls, data: dict) -> 'Params':
        """Create an instance from a JSON file"""
        return cls(**data)

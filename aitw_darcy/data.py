###
#
# Input parameters
#
###

# Label for the simulation
label = '3d_plate'

# Mesh dimensions (m)
xmin = 0
xmax = 5e-2
ymin = 0
ymax = xmax
zmin = 0
zmax = 1e-3
nx = 40             # Number of nodes
ny = 160
nz = 40

# Timestepping
Δt0 = 1e-1                # Initial time step length (s)
Δt_min = 1e-2
Δt_max = 10
t_end = 120               # Simulation time (s)

# Porosity
ϕ = 0.6

# Surface tension (N/m)
γ = 30e-3

# Contact angle
θdeg = 10

# Mean pore radius (m)
r_μ = 20e-6

# Relative standard deviation of pore radius
rstd = 1/5

# Pore radius cutoff
r_cutoff = 2e-6


# Permeability
k_long = 1e-12
k_trans = 1e-16

# Fluid viscosity (Pa s)
μ = 1e-3

# External pressure
u_ext = -1e-6

###
###

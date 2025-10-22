import sys
from runpy import run_path

import numpy as np
import ufl
from basix.ufl import element
from dolfinx import default_scalar_type, fem, io, log, mesh
from dolfinx.fem.petsc import NonlinearProblem
from dolfinx.nls.petsc import NewtonSolver
from matplotlib import pyplot
from mpi4py import MPI
from petsc4py import PETSc
from tqdm import tqdm, trange
from ufl import cos, dot, dx, erf, exp, grad, inner, pi, sqrt

# Inject the variables from the parameter file to the global namespace
try:
    param_file = sys.argv[1]
    params = run_path(param_file)
except:
    error('Please provide the parameter file (.py) in the working directory!')

for k, v in params.items():
    if not k.startswith('_'):
        globals()[k] = v


# Output
solution_file = f'{label}_solution.npy'   # Full solution
output_file = f'{label}_output.bp'        # VTX mesh for visualization

log.set_output_file(f'{label}.log')
log.set_log_level(log.LogLevel.INFO)

# Helper functions

def get_variable_value(v):
    """ Get the nodal values for an expression v """
    q = fem.Function(V)
    q.interpolate(fem.Expression(v, V.element.interpolation_points()))
    return q.x.array

def ramp_up(val=1, t1=0.1):
    """ Linear ramp from zero beginning for t = 0, useful to assist convergence in the beginning """
    return val*ufl.min_value(1, t/t1)

# Generate 3D mesh
domain = mesh.create_box(MPI.COMM_WORLD, [[xmin, ymin, zmin], [xmax, ymax, zmax]], [nx, ny, nz],
                         mesh.CellType.hexahedron)

x = ufl.SpatialCoordinate(domain)

# Current time and time step size
t = fem.Constant(domain, default_scalar_type(0))
Δt = fem.Constant(domain, default_scalar_type(Δt0))

# Create function space
P1 = element('Lagrange', domain.basix_cell(), 1)
V = fem.functionspace(domain, P1)

# Unknown solution (t=n+1)
u = fem.Function(V, name='u')

# Previous time step solution
u_n = fem.Function(V, name='u_n')

# Saturation curve
def S_curve(r):
    σ = rstd*r_μ
    return 1/2*(erf((r-r_μ)/(sqrt(2)*σ)) - erf((r_cutoff-r_μ)/(sqrt(2)*σ)))

# Capillary radius
θ = θdeg*pi/180
r_c = -2*γ*cos(θ)/u

# Saturation
#S = ufl.conditional(u < 0, S_curve(r_c), 1.0)
S = S_curve(r_c)

# Previous saturation
S_n = ufl.replace(S, {u: u_n})

# Flux
k = ufl.as_tensor([[k_long, 0, 0], [0, k_trans, 0], [0, 0, k_trans]])
k_s = 1 + ramp_up(-1+S)
J = -k_s*k/μ*grad(u)

###############
# Weak form
###############

δu = ufl.TestFunction(V)

F = δu * ϕ*(S - S_n) * dx - Δt*dot(grad(δu), J) * dx

# Direchlet boundary conditions

fdim = domain.topology.dim - 1
domain.topology.create_connectivity(fdim, fdim+1)
boundary_facets = mesh.exterior_facet_indices(domain.topology)
boundary_dofs = fem.locate_dofs_topological(V, fdim, boundary_facets)
bc = fem.dirichletbc(default_scalar_type(u_ext), boundary_dofs, V)
bcs = [bc]

# Initial pressure
u0 = -2*γ*cos(θ)/r_cutoff

if __name__ == '__main__':
    # Initial moisture content
    u.interpolate(lambda x: np.full(x[0].shape, u0))
    u.x.array[boundary_dofs] = u_ext

    # Initialize the Newton solver
    problem = NonlinearProblem(F, u, bcs)
    solver = NewtonSolver(MPI.COMM_WORLD, problem)

    # Solution is rather sensitive to tolerance,
    # numerical accuracy seems to be close to 1e-11
    solver.rtol = 0
    solver.atol = 1e-11
    solver.report = True
    solver.error_on_nonconvergence = False
    solver.max_it = 10
    solver.relaxation_parameter = 0.8

    # Use mumps LU solver
    ksp = solver.krylov_solver
    ksp.setType('preonly')
    ksp.getPC().setType('lu')
    ksp.getPC().setFactorSolverType('mumps')

    # Solutions at each time step
    sol = [u.x.array.copy()]

    # Time values
    times = [t.value.item()]

    if domain.comm.rank == 0:
        pbar = tqdm(total=t_end, unit_scale=True)

    S_func = fem.Function(V, name='S')

    with io.VTXWriter(domain.comm, output_file, [u, S_func]) as vtxfile:
        while True:
            u_n.x.array[:] = u.x.array

            # Adaptive timestepping
            while True:
                n, converged = solver.solve(u)
                if converged:
                    solver.relaxation_parameter = 0.8
                    if n <= 3:
                        Δt.value = min(Δt.value * 1.2, Δt_max)
                        log.log(log.LogLevel.INFO, f'Timestep = {Δt.value}')
                    elif n >= 7:
                        Δt.value = Δt.value * 0.5
                        log.log(log.LogLevel.INFO, f'Timestep = {Δt.value}')
                    break
                else:
                    solver.relaxation_parameter = 0.4
                    Δt.value = Δt.value * 0.5
                    if Δt.value < Δt_min:
                        raise Error('Minimum timestep size reached!')
                    u.x.array[:] = u_n.x.array
                    log.log(log.LogLevel.INFO, f'Timestep = {Δt.value}')

            t.value += Δt.value

            # Record values
            sol.append(u.x.array.copy())
            times.append(t.value.item())

            if domain.comm.rank == 0:
                pbar.update(Δt.value)

            # Save output
            S_func.interpolate(fem.Expression(S, V.element.interpolation_points()))
            vtxfile.write(t.value.item())

            # Termination criterion
            if t.value > t_end:
                break

            S_threshold = 0.99
            S_vals = get_variable_value(S)
            finished = np.all(S_vals > S_threshold)
            statuses = domain.comm.gather(finished, root=0)
            can_stop = None
            if domain.comm.rank == 0:
                can_stop = all(statuses)
            if domain.comm.bcast(can_stop, root=0):
                break

    # Numpy array output for single parallel runs
    if domain.comm.size == 1:
        sol = np.array(sol)
        np.save(solution_file, sol)


# The residual and Jacobian for diagnostics
#residual = fem.form(F)
#RD = fem.petsc.create_vector(residual)
#fem.petsc.assemble_vector(RD, residual)
#jacobian = fem.form(ufl.derivative(F, u))

#fem.petsc.apply_lifting(RD, [jacobian], [bcs])

#A = fem.petsc.create_matrix(jacobian)
#fem.petsc.assemble_matrix(A, jacobian)
#A.assemble()
#print(domain.comm.rank, A.getSize())
#m = A.convert('dense').getDenseArray()

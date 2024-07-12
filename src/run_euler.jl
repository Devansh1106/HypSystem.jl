using HypSystem

grid_size = 400     # number of cells

xmin, xmax = 0.0, 1.0       # domain
nvar = 3                    # number of variables
final_time = 0.2
γ = 1.4
disc_x = 0.3                # discontinuity in initial condition in x direction

numflux = "rusanov"

# specify Ul, Ur, in primitive coordinates
# Ul: left of discontinuity; Ur: right side of discontinuity
primitive_l, primitive_r = [1.0, 0.75, 1.0], [0.125, 0.0, 0.1]  # density, velocity, pressure

Ul, Ur = primitive2pde(primitive_l, γ), primitive2pde(primitive_r, γ)
initial_value(x::Float64, eq::Euler) = (x <= disc_x) ? Ul : Ur
boundary_condition = "Periodic"
# initial_value(x::Float64, eq::Euler) = [sin(2.0*pi*x), sin(2.0*pi*x), sin(2.0*pi*x)]
cfl = 0.9

# Parameters printing to screen
println("Solving Riemann problem for Euler equation with following parameters:")
println("Discontinuity of initial data is at ", disc_x)
println("gamma = ", γ)
println("Tf    = ", final_time)
println("Initial Condition L/R of discontinuity ([dens, vel, pres]) is:")
println(primitive_l)
println(primitive_r)

# FVM solver
eq = Euler(γ)
problem = Problem((xmin, xmax), nvar, initial_value, boundary_condition, final_time)
param = create_parameters(cfl, grid_size)
scheme = Scheme(eq, numflux)
grid = make_grid(problem, param)
U = solve(eq, grid, problem, scheme, param)
Ue = compute_exact_soln(eq, grid, final_time, primitive_l, primitive_r, disc_x)
# display(U)
error_cal(eq, grid, U, Ue)

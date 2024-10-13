using HypSystem

grid_size = 800 
xmin, xmax = 0.0, 1.0       # domain
nvar = 3                    # number of variables
fprime = [1.0 2.0 3.0;
          0.0 -2.0 3.0;
          0.0 0.0 3.0]

final_time = 1.0
numflux = "rusanov"
boundary_condition = "Periodic"
initial_value(x::Float64, eq::LinAdv) = [sin(2.0*pi*x),0.5*sin(2.0*pi*x),1.5*sin(2.0*pi*x)]
cfl = 0.9

equation = LinAdv(fprime)
problem = Problem((xmin, xmax), nvar, initial_value, boundary_condition, final_time)
param = create_parameters(cfl, grid_size)
scheme = Scheme(equation, numflux)
grid = make_grid(problem, param)
U = solve(equation, grid, problem, scheme, param)
display(U)

Ue = compute_exact_soln(equation, grid, problem, final_time)
error_cal(equation, grid, U, Ue)
plot_sol(grid, U, Ue, problem)
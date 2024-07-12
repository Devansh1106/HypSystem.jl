module HypSystem

include("grid.jl")
include("hyp_sys1d.jl")

# Linear advection equation
include("EqLinAdv.jl")

# Euler's equation
include("EqEuler.jl")


export compute_exact_soln, error_cal, LinAdv, Euler, primitive2pde, Problem, create_parameters, Scheme, solve, make_grid

end # module HypSystem

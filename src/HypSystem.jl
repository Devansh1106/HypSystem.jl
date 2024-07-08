module HypSystem

include("grid.jl")
include("hyp_sys1d.jl")
include("EqLinAdv.jl")

export LinAdv, Problem, create_parameters, Scheme, solve

end # module HypSystem

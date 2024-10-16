struct CartesianGrid
    domain::Tuple{Float64, Float64} # xmin, xmax
    nx::Int64                   # nx - number of points
    xc::Vector{Float64}         # cell centers
    xf::Vector{Float64}         # cell faces
    dx::Vector{Float64}         # cell sizes
    nvar::Int64
end

# Uniform cartesian grid 
function make_grid(problem, param)
    nvar = problem.nvar
    xmin, xmax = problem.domain 
    nx = param.grid_size
    println("Making Uniform Grid of interval [$xmin, $xmax]")
    dx = (xmax - xmin)/nx
    xc = LinRange(xmin+0.5*dx, xmax-0.5*dx, nx)
    println("Grid with number of points = $nx")
    println("(xmin, xmax) = ($xmin, $xmax)")
    xf = LinRange(xmin, xmax, nx+1)
    dx1 = fill(0.0, nx)
    dx1 .= dx
    return CartesianGrid((xmin, xmax), nx, xc, xf, dx1, nvar)
end

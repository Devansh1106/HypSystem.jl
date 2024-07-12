using LinearAlgebra
using OffsetArrays
using DelimitedFiles
using Plots

struct Problem{F <: Function}
    domain::Tuple{Float64, Float64}
    nvar::Int64
    initial_value::F
    # boundary_value::F
    boundary_condition::Any
    final_time::Float64
end

struct Scheme
    equation
    numflux::String
end

struct Parameters
    cfl::Float64
    grid_size::Int64
end

function create_parameters(cfl::Float64, grid_size::Int64)
    @assert (cfl >= 0.0) "cfl should be >= 0.0"
    return Parameters(cfl, grid_size)
end

gArray(nvar::Int64, nx::Int64) = OffsetArray(zeros(nvar, nx+2),
                                             OffsetArrays.Origin(1,0))

function adjust_time_step(problem::Problem, param::Parameters,
                          dt::Float64, t::Float64)
    final_time = problem.final_time
    if t + dt > final_time
        dt = final_time - t
    end

    return dt    
end

# TODO: for non-uniform grid there will be loop and dt = min(...);lam=max(...)
function compute_lam_dt(U::gArray, param::Parameters, grid::CartesianGrid, scheme::Scheme, eq) where gArray
    nx = grid.nx
    dx = grid.dx
    lam = 0.0
    dt = 1.0
    for i in 1:nx
        ua = U[:, i]
        if typeof(eq) == LinAdv
            lam0 = maximum(abs.(eigvals(eq.fprime)))
        else
            ρ, u, E = ua[1], ua[2]/ua[1], ua[3]    # density, velocity, energy
            p = (eq.γ - 1.0) * (E - 0.5*ρ*u^2)     # pressure
            c = sqrt(eq.γ*p/ρ)                     # sound speed
            # @show eq.γ*p/ρ
            lam0 = abs(u)+c
        end
        lam = max(lam, lam0)
        dt = min(dt, dx[i]/lam0)
     end

     dt = param.cfl * dt   # dx is a vector; it is useful for non-uniform grid version
     lam = param.cfl * lam
    return lam, dt
end

function set_initial_value!(eq, grid::CartesianGrid, U::gArray, problem::Problem) where gArray
    nx = grid.nx
    xc = grid.xc
    initial_value = problem.initial_value
    for i in 1:nx
        @views U[:,i] .= initial_value(xc[i], eq)
    end
end

function update_ghost!(grid::CartesianGrid, U::gArray, problem::Problem) where gArray
    nx = grid.nx
    if problem.boundary_condition == "Periodic"
        U[:, 0] .= U[:, nx]
        U[:, nx+1] .= U[:, 1]
    end
end

function compute_residual!(equation, lam::Float64, grid::CartesianGrid, U::gArray,
                            scheme::Scheme, res::gArray, problem::Problem) where gArray
    Uf = zeros(problem.nvar)
    nx = grid.nx
    dx = OffsetArray(zeros(nx+2), OffsetArrays.Origin(0))
    @. dx[1:nx] = grid.dx[1:nx]
    res[:,:] .= 0.0     # Every time step we need to make res = 0.0, otherwise it will add up all
    # dx[0] = grid.dx[nx]	# no need of this as res[0] and res[nx+1] are of no use.
    # dx[nx+1] = grid.dx[1]
    for i in 1:nx+1
        @views Ul, Ur = U[:, i-1], U[:, i]
        if scheme.numflux == "rusanov"
           rusanov!(equation, lam, Ul, Ur, Uf)
        end
        @views res[:, i-1] += Uf/ dx[i-1]           
        @views res[:, i] -= Uf/ dx[i]
    end
end

function solve(equation, grid::CartesianGrid, problem::Problem, scheme::Scheme, param::Parameters)
    nvar = problem.nvar
    Tf = problem.final_time
    nx = grid.nx
    dx = grid.dx
    xf = grid.xf
    # Allocating variables
    U = gArray(nvar, nx)
    res = gArray(nvar, nx) # dU/dt + res(U) = 0
    # display(U)
    set_initial_value!(equation, grid, U, problem)
    it, t = 0, 0.0
    while t < Tf 
        lam, dt = compute_lam_dt(U, param, grid, scheme, equation)
        dt = adjust_time_step(problem, param, dt, t)
        update_ghost!(grid, U, problem)
        compute_residual!(equation, lam, grid, U, scheme, res, problem)
        @. U -= dt*res
        t += dt; it += 1
    end
    return U
   # if typeof(equation) == LinAdv
   #     compute_exact_soln!(grid, equation, problem, t, Ue)
   # else
   #     println("calculating the exact solution from Toro's solver")
   #     println("Solving exactly for final time")
   #     p_l_s, p_r_s = array2string(primitive_l), array2string(primitive_r)
   #     # run(`python3 ./ToroExact/toro_exact.py -p user -l $p_l_s -r $p_r_s -x $disc_x -t $final_time`)
   #     print("python3 ./ToroExact/toro_exact.py -p user -l $p_l_s -r $p_r_s -x $disc_x -t")
   # end
    # error_cal!(grid, problem, U, Ue, error) 
    # display(error)
    # plot_sol(grid, U, Ue, problem)
end

function plot_sol(grid::CartesianGrid, U::gArray, Ue::gArray, problem::Problem) where gArray
    nvar = problem.nvar
    nx = grid.nx
    xc = grid.xc
    for i in 1:nvar
        @views y_exact, y_num = Ue[i, 1:nx], U[i, 1:nx]
        plot(xc, y_exact, label="Exact Solution", linestyle=:solid, linewidth=2, dpi=150)
        plot!(xc, y_num, label="Numerical Solution", xlabel="Domain", ylabel="Solution values",
              title="Solution Plot $nvar", linewidth=2, linestyle=:dot, linecolor="black", 
              dpi=150)
        savefig("fig/hypsys1D_$i.png")
    end
    println("Figures are saved in folder `fig`.")
end

# function error_cal!(grid::CartesianGrid, problem::Problem, U::gArray, Ue::gArray, error::Vector{Float64}) where gArray
#     nx = grid.nx
#     nvar = problem.nvar
#     for i in 1:nvar
# 	@views error[i] = sum(abs.(U[i, 1:nx] - Ue[i, 1:nx]))
# 	error[i] = error[i]/nx
#     end
# end

# converts array to string; used in Toro's solver for Euler's equation
function array2string(arr)
   arr_string = "["
   n = size(arr)[1]
   for i=1:n-1
      arr_string = arr_string * string(arr[i]) * ","
   end
   arr_string = arr_string * string(arr[end]) * "]"
end

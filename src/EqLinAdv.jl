using LinearAlgebra
using DelimitedFiles

struct LinAdv
    fprime::Array{Float64, 2}
end

function flux(fprime::Array{Float64, 2}, U::SubArray{Float64, 1})
    F = fprime * U
    return F
end

function rusanov!(equation::LinAdv, lam::Float64, Ul::SubArray{Float64, 1},
                  Ur::SubArray{Float64, 1}, Uf::Vector{Float64})
    lam = maximum(abs.(eigvals(equation.fprime)))
    fprime = equation.fprime
    Fl, Fr = flux(fprime, Ul), flux(fprime, Ur)
    Uf .= 0.5*(Fl + Fr) - 0.5*lam*(Ur - Ul)
end

function compute_exact_soln(equation::LinAdv, grid::CartesianGrid, problem::Problem,
                             t::Float64)
    nx = grid.nx
    xc = grid.xc
    nvar = problem.nvar
    Ue = gArray(nvar, nx)
    initial_value = problem.initial_value
    eigen_vals = eigvals(equation.fprime)
    eigen_vecs = eigvecs(equation.fprime)
    for j in 1:nx
        for i in 1:nvar
            Ue[i,j] = (inv(eigen_vecs) * initial_value(xc[j] - eigen_vals[i] * t, equation))[i]
        end
    end
    for i in 1:nx
        @views Ue[:,i] .= eigen_vecs * Ue[:, i]
    end
    return Ue
end

function error_cal(eq::LinAdv, grid::CartesianGrid, U::gArray, Ue::gArray) where gArray
    nx = grid.nx
    nvar = grid.nvar
    error = fill(0.0, nvar)
    for i in 1:nvar
        @views error[i] = sum(abs.(U[i, 1:nx] - Ue[i, 1:nx]))
        error[i] = error[i]/nx
    end
    display(error)
end


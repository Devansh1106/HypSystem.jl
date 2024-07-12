using LinearAlgebra
using DelimitedFiles

struct Euler
    γ::Float64
end

function flux(U::SubArray{Float64, 1}, eq::Euler)
    ρ = U[1]        # density
    u = U[2]/U[1]   # velocity
    E = U[3]        # energy
    p = (eq.γ - 1.0) * (E - 0.5 * ρ * u^2) # pressure
    F = [U[2], p + ρ * u^2, (E+p) * u] # flux
    return F
end

# Numerical flux
function rusanov!(eq::Euler, lam::Float64, Ul::SubArray{Float64, 1}, Ur::SubArray{Float64, 1},
                  Uf::Vector{Float64})
   γ = eq.γ
   ρl, ul, pl = pde2primitive(Ul, γ)
   ρr, ur, pr = pde2primitive(Ur, γ)
   cl, cr = sqrt(γ*pl/ρl), sqrt(γ*pr/ρr)                   # sound speed
   λ = maximum(abs.([ul, ul-cl, ul+cl, ur, ur-cr, ur+cr])) # local wave speed
   Fl, Fr = flux(Ul, eq), flux(Ur, eq)
   Uf  .= 0.5*(Fl+Fr) - 0.5*λ*(Ur - Ul)
end

# converting primitive variables to PDE variables
function primitive2pde(prim::Vector{Float64}, γ::Float64)
    U = [prim[1], prim[1]*prim[2], prim[3]/(γ-1.0)/prim[1] + (prim[1]*prim[2]^2) /2.0]
   #     ρ     ,   ρ*u          , p/(γ-1.0)/ρ + ρ*u^2 /2.0 = E
   return U
end

# converting PDE variables to primitive variables
function pde2primitive(U::SubArray{Float64, 1}, γ::Float64)
    primitives = [U[1], U[2]/U[1], U[1]*(γ-1.0)*(U[3]-U[2]^2/2.0)]
   #            [ρ ,      u          , p]
   return primitives
end

function error_cal(eq::Euler, grid::CartesianGrid, U::gArray, Ue::gArray) where gArray
    nx = grid.nx
    nvar = grid.nvar
    error = fill(0.0, nvar)
    display(U)
    display(Ue)
    for i in 1:nvar
        @views error[i] = sum(abs.(U[i, 1:nx] - Ue[i, 1:nx]))
        error[i] = error[i]/nx
    end
    display(error)
end

function compute_exact_soln(eq::Euler, grid::CartesianGrid, t::Float64, primitive_l::Vector{Float64}, primitive_r::Vector{Float64}, disc_x::Float64)
    nx = grid.nx
    nvar = grid.nvar
    Ue = gArray(nvar, nx)
    println("calculating the exact solution from Toro's solver")
    println("Solving exactly for final time")
    p_l_s, p_r_s = array2string(primitive_l), array2string(primitive_r)
    run(`python3 ./ToroExact/toro_exact.py -p user -n $nx -l $p_l_s -r $p_r_s -x $disc_x -t $t`)
    # print("python3 ./ToroExact/toro_exact.py -p user -l $p_l_s -r $p_r_s -x $disc_x -t")
    exact_data = readdlm("./ToroExact/output/toro_user_exact.dat")
    # display(exact_data[10:grid_size, 2])
    for i in 1:nx
        Ue[1, i] = exact_data[9+i, 2]
        Ue[2, i] = exact_data[9+i, 4]
        Ue[3, i] = exact_data[9+i, 3]
    end
    return Ue
end

# Libraries
using LinearAlgebra
using Distributions
using Interpolations
using Optim
using Plots
using StatsBase
using QuantEcon
using LaTeXStrings
import Base.@kwdef

plot_font = "Computer Modern"
default(fontfamily=plot_font,
        linewidth=2, framestyle=:box, label=nothing, grid=false)
scalefontsizes(1)

# Datastruct for the model 
@kwdef struct Model_Params

    # Asset grid paramers
    agrid::Int16 = 50
    amin::Float32 = 0.01
    amax::Float32 = 5.0
    A = LinRange(amin, amax, agrid)

    # Idiosyncratic Shock
    egrid::Int16 = 30
    ρ::Float32 = 0.7
    nstd::Int16 = 2
    σ_ρ::Float32 = 0.5
    E = QuantEcon.tauchen(egrid, ρ, σ_ρ, 0, nstd).state_values
    Π = QuantEcon.tauchen(egrid, ρ, σ_ρ, 0, nstd).p

    # Occupations
    ogrid::Int16 = 3
    σ_ϕ::Int16 = 1
    O::LinRange{Int64, Int64} = LinRange(1, ogrid, ogrid)

    # Utility and Preferences
    σ::Int16 = 2
    β::Float32 = 0.96

    # Productivity and types
    Γ = [1 0.32; 
         0.46 1]

    m1::Float32 = 0.7
    m2::Float32 = 0.3

    # Firm
    α::Float32 = 0.36
    ξ::Float32 = 0.36
    price_1::Float32 = 1.0
    price_2::Float32 = 1.0

    # Prices
    r::Float32 = 0.03
    u::Float32 = 0.1

    # Tolerance
    ε = 1e-04
end

function ufun(c, σ)
    if σ == 1 
        return log(c) 
    else 
        return (c^(1 - σ))/(1 - σ) 
    end
end

function find_closest(arr, target)
    # Calculate the absolute differences between each element and the target
    differences = abs.(arr .- target)
    
    # Find the indices that would sort the differences
    sorted_indices = sortperm(differences)
    
    # Take the two closest values based on the sorted indices
    closest_values = arr[sorted_indices[1:2]]
    
    return closest_values
end

function Φ(params::Model_Params, AP)
    (;agrid, egrid, A, E, Π) = params
    next::Array{Float64, 4} = ones(agrid, egrid, agrid, egrid)

    Threads.@threads for idx1 in eachindex(A) # a
                    for idx5 in eachindex(E) # ε
                        Threads.@threads for idx6 in eachindex(A) # a'
                                        for idx10 in eachindex(E) # ε'
                                            closest = find_closest(A, AP[idx1, idx5])
                                            @inbounds next[idx1, idx5, idx6, idx10] = (A[idx6] == closest[1])*Π[idx5, idx10]
                                        end
                                    end
                                end
                            end

    Pmat = reshape(next, (agrid*egrid, agrid*egrid))
    ϕ = QuantEcon.stationary_distributions(QuantEcon.MarkovChain(reverse(Pmat)))
    return reshape(ϕ[1], (agrid, egrid))
end

function wealth_distribution(params::Model_Params, λ)
    (; agrid) = params
    return reshape(sum(λ, dims=(2)), agrid)
end

function ∫(EV, Ap, Ep, CurrE, Π)
    logsum = EV(Ap, Ep)'*Π[CurrE, :]
    return logsum
end

function T(params::Model_Params, r,  w1, w2, u, EV, Γ)

    (; agrid, egrid, ogrid, β, σ_ϕ, σ, A, E, Π) = params

    # Istantiate arrays
    v::Array{Float64, 3} = zeros(agrid, egrid, ogrid)
    c::Array{Float64, 3} = zeros(agrid, egrid, ogrid)

    # For each point in the domain
    Threads.@threads for idx1 in eachindex(A) # ∀ a ∈ A
        @inbounds for idx5 in eachindex(E)  # ∀ ε ∈ E
            
            # Choice 1
            v1 = c -> -(ufun(c, σ) + β*∫(EV, ((1+r)*A[idx1] + Γ[1]*exp(E[idx5])*w1 - c), E, idx5, Π))

            # Choice 2
            v2 = c -> -(ufun(c, σ) + β*∫(EV, ((1+r)*A[idx1] + Γ[2]*exp(E[idx5])*w2 - c), E, idx5, Π))

            # Choice 3
            v3 = c -> -(ufun(c, σ) +  β*∫(EV, ((1+r)*A[idx1] + u - c), E, idx5, Π))
            
            # Find solution to optimization
            res1 = optimize(v1, 0, (1+r)*A[idx1] + Γ[1]*exp(E[idx5])*w1)
            res2 = optimize(v2, 0, (1+r)*A[idx1] + Γ[2]*exp(E[idx5])*w2)
            res3 = optimize(v3, 0, (1+r)*A[idx1] + u)

            # Instantiate the points
            v[idx1, idx5, 1] = -Optim.minimum(res1)
            v[idx1, idx5, 2] = -Optim.minimum(res2)
            v[idx1, idx5, 3] = -Optim.minimum(res3)

            c[idx1, idx5, 1] = Optim.minimizer(res1)[1]
            c[idx1, idx5, 2] = Optim.minimizer(res2)[1]
            c[idx1, idx5, 3] = Optim.minimizer(res3)[1]
        end
    end

    # Numerical stability constant
    M = maximum(v)

    # Compute the logsum 
    EV1 = M .+ σ_ϕ.*log.(exp.((v[:,:,1] .- M)./σ_ϕ) .+ exp.((v[:,:,2] .- M)./σ_ϕ) .+ exp.((v[:,:,3] .- M)./σ_ϕ))

    # Interpolation step
    EV1 = linear_interpolation((A, E), EV1, extrapolation_bc=Line());

    return EV1, v, c
end

function VFI(params::Model_Params, type, r, w1, w2, u, print_res)

    (; agrid, A, egrid, E, ogrid, ε, Γ) = params

    if type == 1
        Γ = Γ[1,:]
    else
        Γ = Γ[2,:]
    end

    # Initial Guess on EV
    EV = linear_interpolation((A, E), zeros(agrid, egrid), extrapolation_bc=Line());

    # Instantiate v and c
    v::Array{Float64, 3} = zeros(agrid, egrid, ogrid)
    c::Array{Float64, 3} = zeros(agrid, egrid, ogrid)

    if print_res == true
        println("Starting Value Function Iteration")
    end
    err = Inf;
    # Solution loop
    while err > ε
        EV1, v, c = T(params, r,  w1, w2, u, EV, Γ)
        err = StatsBase.Linfdist(EV1(A, E), EV(A, E))
        if print_res == true
            println(err)
        end
        EV = EV1
    end

    return EV, v, c
end

function CCP(params::Model_Params, v)
    (;σ_ϕ) = params
    M = maximum(v)
    p1 = exp.((v[:,:,1] .- M)./σ_ϕ)./(exp.((v[:,:,1] .- M)./σ_ϕ) + exp.( (v[:,:,2] .- M)./σ_ϕ) + exp.((v[:,:,3] .- M)./σ_ϕ))
    p2 = exp.((v[:,:,2] .- M)./σ_ϕ)./(exp.((v[:,:,1] .- M)./σ_ϕ) .+ exp.((v[:,:,2] .- M)./σ_ϕ) .+ exp.((v[:,:,3] .- M)/σ_ϕ))
    pu = 1 .- p1 .- p2
    return p1, p2, pu
end

function ap(params::Model_Params, type, r, w1, w2, u, p1, p2, pu, c)

    (;agrid, egrid, ogrid, A, E, Γ) = params

    if type == 1
        Γ = Γ[1,:]
    else
        Γ = Γ[2,:]
    end

    ap::Array{Float64, 3} = zeros(agrid, egrid, ogrid)

    Threads.@threads for idx1 in eachindex(A)
        @inbounds for idx5 in eachindex(E)
            ap[idx1, idx5, 1] = (1 + r)*A[idx1] + Γ[1]*exp(E[idx5])*w1 - c[idx1, idx5, 1]
            ap[idx1, idx5, 2] = (1 + r)*A[idx1] + Γ[2]*exp(E[idx5])*w2 - c[idx1, idx5, 2]
            ap[idx1, idx5, 3] = (1 + r)*A[idx1] + u - c[idx1, idx5, 3]
        end
    end

    return ap[:,:,1].*p1 + ap[:,:,2].*p2 + ap[:,:,3].*pu;
end

# Solve the Model
function solve_infinite_horizon(params, type, r, w1, w2, u, print_res)
    # Value Function Iteration
    EV, v, c = VFI(params, type, r, w1, w2, u, print_res)

    # Compute choice probabilities
    p1, p2, pu = CCP(params, v)

    # Reconstruct ap
    AP = ap(params, type, r, w1, w2, u, p1, p2, pu, c)

    C = c[:,:,1].*p1 + c[:,:,2].*p2 + c[:,:,3].*pu;
    
    λ = Φ(params, AP)

    WD = wealth_distribution(params, λ);

    return EV, p1, p2, pu, v, C, AP, λ, WD
end

# General equilibrium algorithms
function wages(params::Model_Params, L1, L2)
    (; α, ξ, price_1, price_2) = params
    w1 = α*price_1*L1^(α-1)
    w2 = ξ*price_2*L2^(ξ-1)
    return w1, w2
end

function find_equilibrium(params::Model_Params, L1_init, L2_init, b_weights)
    (; Γ, r, u) = params
    ε = 1e-02
    it = 0
    it_max = 10

    # Inital guess of capital and labor
    L1 = L1_init
    L2 = L2_init

    w1 = 0 
    w2 = 0

    while it < it_max
        w1, w2 = wages(params, L1, L2)
        println("Current prices, r=$r, w1=$w1, w2=$w2")

        # Compute the model at that eq_prices
        EV_1, p1_1, p2_1, pu_1, v_1, C_1, AP_1, λ_1, WD_1 = solve_infinite_horizon(params, 1, r, w1, w2, u, false) # Solve the Model for type 1
        EV_2, p1_2, p2_2, pu_2, v_2, C_2, AP_2, λ_2, WD_2 = solve_infinite_horizon(params, 2, r, w1, w2, u, false) # Solve the model for type 2

        # Compute aggregate capital and labor supply
        L1_new = effective_labor_supply(params, p1_1, p1_2, λ_1, λ_2, [Γ[1,1], Γ[2,1]])
        L2_new = effective_labor_supply(params, p2_1, p2_2, λ_1, λ_2, [Γ[1,2], Γ[2,2]])

        println(abs(L1_new - L1))
        println(abs(L2_new - L2))

        if (abs(L1_new - L1) < ε) & (abs(L2_new - L2) < ε)
            println("Approximate equilibrium found at r=$r, w1 = $w1, w2 = $w2")
            break
        else
            L1 = b_weights*L1_new + (1 - b_weights)*L1
            L2 = b_weights*L2_new + (1 - b_weights)*L2 
            it = it + 1
        end
    end

    return w1, w2, r, u
end

function employment_rate(params::Model_Params, p1, p2, λ1, λ2)
    (;m1, m2) = params
    E11 = p1.*λ1
    E12 = p2.*λ2
    return sum(m1*E11 + m2*E12), sum(E11), sum(E12)
end

function effective_labor_supply(params::Model_Params, p1, p2, λ1, λ2, Γ)
    (;m1, m2) = params
    E11 = p1.*λ1
    E12 = p2.*λ2
    return sum(m1*E11*Γ[1] + m2*E12*Γ[2])
end

function 𝒞(params::Model_Params, c1, c2, λ1, λ2)
    (;m1, m2) = params
    C1 = c1.*λ1
    C2 = c2.*λ2
    return sum(m1*C1 + m2*C2)
end

function 𝒦(params::Model_Params, AP1, AP2, λ1, λ2)
    (;m1, m2) = params
    S1 = AP1.*λ1
    S2 = AP2.*λ2
    return sum(m1*S1 + m2*S2)
end

function sectoral_shock(params::Model_Params, w1, w2lb, w2ub, n_grid)
    
    (;r, u, m1, m2, agrid, Γ) = params

    W2 = LinRange(w2lb, w2ub, n_grid)
    WDIST = zeros(n_grid, agrid)
    summary_mat = zeros(18, n_grid)

    Threads.@threads for i in eachindex(W2) # For every wage 2, in parallel

        # Compute the model at that eq_prices
        EV_1, p1_1, p2_1, pu_1, v_1, C_1, AP_1, λ_1, WD_1 = solve_infinite_horizon(params, 1, r, w1, W2[i], u, false) # Solve the Model for type 1
        EV_2, p1_2, p2_2, pu_2, v_2, C_2, AP_2, λ_2, WD_2 = solve_infinite_horizon(params, 2, r, w1, W2[i], u, false) # Solve the model for type 2

        # Compute aggregate measures
        E1, E11, E21 = employment_rate(params, p1_1, p1_2, λ_1, λ_2)
        E2, E12, E22 = employment_rate(params, p2_1, p2_2, λ_1, λ_2)
        E = E1 + E2
        U, U1, U2 =  employment_rate(params, pu_1, pu_2, λ_1, λ_2)

        L1 = effective_labor_supply(params, p1_1, p1_2, λ_1, λ_2, [Γ[1,1], Γ[2,1]])
        L2 = effective_labor_supply(params, p2_1, p2_2, λ_1, λ_2, [Γ[1,2], Γ[2,2]])
        L = L1 + L2
        
        mean_prod_1 = L1/E1
        mean_prod_2 = L2/E2
        prod = L/E

        C = 𝒞(params, C_1, C_2, λ_1, λ_2) 
        K = 𝒦(params, AP_1, AP_2, λ_1, λ_2)

        summary_mat[1, i] = E
        summary_mat[2, i] = U
        summary_mat[3, i] = E1
        summary_mat[4, i] = E2
        summary_mat[5, i] = E11
        summary_mat[6, i] = E21
        summary_mat[7, i] = E12
        summary_mat[8, i] = E22
        summary_mat[9, i] = U1
        summary_mat[10, i] = U2

        summary_mat[11, i] = L1
        summary_mat[12, i] = L2
        summary_mat[13, i] = L

        summary_mat[14, i] = mean_prod_1
        summary_mat[15, i] = mean_prod_2
        summary_mat[16, i] = prod

        summary_mat[17, i] = C
        summary_mat[18, i] = K

        WDIST[i,:] = m1*WD_1 + m2*WD_2
    end

    return summary_mat, W2, WDIST
end


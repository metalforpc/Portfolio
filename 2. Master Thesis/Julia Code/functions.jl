using LinearAlgebra
using Distributions
using Interpolations
using Optim
using Plots
using StatsBase
using LaTeXStrings
using QuantEcon
using IJulia
using Inequality

plot_font = "Computer Modern"
default(fontfamily=plot_font,
        linewidth=2, framestyle=:box, label=nothing, grid=false)
scalefontsizes(1)

function params()
    # Asset grid
    agrid = 50
    amin = 0.01
    amax = 5
    A = LinRange(amin, amax, agrid)

    # Idiosyncratic Shock
    egrid = 30
    ρ = 0.7
    nstd = 2
    σ_ρ = 0.5
    E = QuantEcon.tauchen(egrid, ρ, σ_ρ, 0, nstd).state_values
    Π = QuantEcon.tauchen(egrid, ρ, σ_ρ, 0, nstd).p

    # Occupations
    ogrid = 3
    σ_ϕ = 1
    O::LinRange{Int64, Int64} = LinRange(1, ogrid, ogrid)

    # Utility and Preferences
    σ = 2
    β = 0.96

    return agrid, A, egrid, E, Π, o, O, σ, σ_ϕ, β
end

function find_closest_numbers(arr, target)
    # Calculate the absolute differences between each element and the target
    differences = abs.(arr .- target)
    
    # Find the indices that would sort the differences
    sorted_indices = sortperm(differences)
    
    # Take the two closest values based on the sorted indices
    closest_values = arr[sorted_indices[1:2]]
    
    return closest_values
end

function stationary_dist(A, E, Π, AP)
    next::Array{Float64, 4} = ones(length(A), length(E), length(A), length(E))

    Threads.@threads for idx1 in eachindex(A) # a
                    for idx5 in eachindex(E) # ε
                         for idx6 in eachindex(A) # a'
                                        for idx10 in eachindex(E) # ε'
                                            closest = find_closest_numbers(A, AP[idx1, idx5])
                                            @inbounds next[idx1, idx5, idx6, idx10] = (A[idx6] == closest[1])*Π[idx5, idx10]
                                        end
                                    end
                                end
                            end

    #println("Transition Matrix computed")

    # Once we computed our matrix

    Pmat = reshape(next, (length(A)*length(E), length(A)*length(E)))
    stationary_distribution = QuantEcon.stationary_distributions(QuantEcon.MarkovChain(reverse(Pmat)))
    n_distr = length(stationary_distribution)
    #println("There are $n_distr stationary distributions, picking the first")
    
    return reshape(stationary_distribution[1], (length(A), length(E)))
end

function wealth_distribution(λ, A, E)
    return reshape(sum(λ, dims=(2)),length(A))
end

function ufun(c, σ)
    if σ != 1
        return (c.^(1 - σ) - 1)./(1 - σ);
    else 
        return log.(c)
    end
end

function intEV(EV, Ap, Ep, CurrE, Π)
    logsum = EV(Ap, Ep)'*Π[CurrE, :]
    return logsum
end

function T(agrid, egrid, o, r, β, σ_ϕ, σ, w1, w2, u, A, E, O, Π, EV, Γ)

    v::Array{Float64, 3} = zeros(agrid, egrid, o)
    c::Array{Float64, 3} = zeros(agrid, egrid, o)

    # For each point in the domain
    Threads.@threads for idx1 in eachindex(A) # ∀ a ∈ A
                    for idx5 in eachindex(E)  # ∀ ε ∈ E
                        
                        # Choice 1
                        v1 = c -> -(ufun(c, σ) + β*intEV(EV, ((1+r)*A[idx1] + Γ[1]*exp(E[idx5])*w1 - c), E, idx5, Π))
            
                        # Choice 2
                        v2 = c -> -(ufun(c, σ) + β*intEV(EV, ((1+r)*A[idx1] + Γ[2]*exp(E[idx5])*w2 - c), E, idx5, Π))

                        # Choice 3
                        v3 = c -> -(ufun(c, σ) +  β*intEV(EV, ((1+r)*A[idx1] + u - c), E, idx5, Π))
                        
                        # Find solution to optimization
                        res1 = optimize(v1, 0, (1+r)*A[idx1] + Γ[1]*exp(E[idx5])*w1)
                        res2 = optimize(v2, 0, (1+r)*A[idx1] + Γ[2]*exp(E[idx5])*w2)
                        res3 = optimize(v3, 0, (1+r)*A[idx1] + u)

                        # Instantiate the points
                        @inbounds v[idx1, idx5, 1] = -Optim.minimum(res1)
                        @inbounds v[idx1, idx5, 2] = -Optim.minimum(res2)
                        @inbounds v[idx1, idx5, 3] = -Optim.minimum(res3)

                        @inbounds c[idx1, idx5, 1] = Optim.minimizer(res1)[1]
                        @inbounds c[idx1, idx5, 2] = Optim.minimizer(res2)[1]
                        @inbounds c[idx1, idx5, 3] = Optim.minimizer(res3)[1]
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

# Solve the Model
function solver(r, w1, w2, u, Γ)

    # Parameters
    agrid, A, egrid, E, Π, o, O, σ, σ_ϕ, β = params

    # Initial Guess on EV
    EV = linear_interpolation((A, E), zeros(agrid, egrid), extrapolation_bc=Line());

    # Instantiate v and c
    v::Array{Float64, 3} = zeros(agrid, egrid, o)
    c::Array{Float64, 3} = zeros(agrid, egrid, o)

    # Tolerance 
    ε = 1e-04;
    err = Inf;

    println("Starting Value Function Iteration")

    # Solution loop
    while err > ε

        EV1, v, c = T(agrid, egrid, o, r, β, σ_ϕ, σ, w1, w2, u, A, E, O, Π, EV, Γ)
    
        err = StatsBase.Linfdist(EV1(A, E), EV(A, E))

        EV = EV1

    end

    # Compute choice probabilities
    M = maximum(v)
    p1 = exp.((v[:,:,1] .- M)./σ_ϕ)./(exp.((v[:,:,1] .- M)./σ_ϕ) + exp.( (v[:,:,2] .- M)./σ_ϕ) + exp.((v[:,:,3] .- M)./σ_ϕ))
    p2 = exp.((v[:,:,2] .- M)./σ_ϕ)./(exp.((v[:,:,1] .- M)./σ_ϕ) .+ exp.((v[:,:,2] .- M)./σ_ϕ) .+ exp.((v[:,:,3] .- M)/σ_ϕ))
    pu = 1 .- p1 .- p2

    # Reconstruct ap
    ap::Array{Float64, 3} = zeros(agrid, egrid, o)
    for (idx1, a) in enumerate(A)
        for (idx5, e) in enumerate(E)
            ap[idx1, idx5, 1] = (1 + 0.03)*a + Γ[1]*exp(E[idx5])*w1 - c[idx1, idx5, 1]
            ap[idx1, idx5, 2] = (1 + 0.03)*a + Γ[2]*exp(E[idx5])*w2 - c[idx1, idx5, 2]
            ap[idx1, idx5, 3] = (1 + 0.03)*a + u - c[idx1, idx5, 3]
        end
    end

    C = c[:,:,1].*p1 + c[:,:,2].*p2 + c[:,:,3].*pu;
    AP = ap[:,:,1].*p1 + ap[:,:,2].*p2 + ap[:,:,3].*pu;

    λ = stationary_dist(A, E, Π, AP)

    WD = wealth_distribution(λ, A, E);

    return EV, p1, p2, pu, v, A, E, C, AP, λ, WD
end

function solve_backward(T, init_L1, init_L2, EV0, EVT, agrid, r, p1, p2)

    # Instantiate the array for the value functions
    EV_tab = zeros(1:T)
    EV_tab[1] = EV0
    EV_tab[T] = EVT

    # Arrays for aggregate consumption and prices
    C = zeros(1:T)
    c = zeros(1:T)
    v = zeros(1:T)
    W1 = zeros(1:T)
    W2 = zeros(1:T)

    # For each period in time, solve backward the model
    for t in reverse(3:T)

        # Compute prices
        w1 = α*p1*init_L1[t]^(α-1)
        w2 = ξ*p2*init_L2[t]^(ξ-1)
        W1[t] = w1
        W2[t] = w2

        # Solve backward recursion
        EV_tab[t-1], v[t-1], c[t-1] = T(agrid, egrid, o, r, β, σ_ϕ, σ, w1, w2, u, A, E, O, Π, EV_tab[t], Γ)

        # Once the backward recursion is solved compute CPP
        M = maximum(v)
        p1 = exp.((v[:,:,1] .- M)./σ_ϕ)./(exp.((v[:,:,1] .- M)./σ_ϕ) + exp.( (v[:,:,2] .- M)./σ_ϕ) + exp.((v[:,:,3] .- M)./σ_ϕ))
        p2 = exp.((v[:,:,2] .- M)./σ_ϕ)./(exp.((v[:,:,1] .- M)./σ_ϕ) .+ exp.((v[:,:,2] .- M)./σ_ϕ) .+ exp.((v[:,:,3] .- M)/σ_ϕ))
        pu = 1 .- p1 .- p2
    end

    
end

function wages(α, ξ, p1, p2, L1, L2)
    w1 = α*p1*L1^(α-1)
    w2 = ξ*p2*L2^(ξ-1)
    return w1, w2
end

########################################################################################## GE ROUTINE
function eq_prices(α, ξ, r, u, m1, m2, L1_init, L2_init, b_weights, p1, p2)

    # Parameters
    Γ_1 = [1 0.32]
    Γ_2 = [0.45 1]

    ε = 1e-02
    it = 0
    it_max = 10

    # Inital guess of capital and labor
    L1 = L1_init
    L2 = L2_init

    w1 = 0 
    w2 = 0

    while it < it_max
        w1, w2 = wages(α, ξ, p1, p2, L1, L2)
        println("Current prices, r=$r, w1=$w1, w2=$w2")

        # Compute the model at that eq_prices
        EV_1, p1_1, p2_1, pu_1, v_1, A, E, C_1, AP_1, λ_1, WD_1 = solver(r, w1, w2, u, Γ_1) # Solve the Model for type 1
        EV_2, p1_2, p2_2, pu_2, v_2, A, E, C_2, AP_2, λ_2, WD_2 = solver(r, w1, w2, u, Γ_2) # Solve the model for type 2

        # Compute aggregate capital and labor supply
        L1_new = effective_labor_supply(p1_1, p1_2, m1, m2, λ_1, λ_2, A, E, [Γ_1[1], Γ_2[1]])
        L2_new = effective_labor_supply(p2_1, p2_2, m1, m2, λ_1, λ_2, A, E, [Γ_1[2], Γ_2[2]])

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

function employment_rate(p1, p2, m1, m2, λ1, λ2, A, E)
    E11 = p1.*λ1
    E12 = p2.*λ2
    return sum(m1*E11 + m2*E12)
end

function effective_labor_supply(p1, p2, m1, m2, λ1, λ2, A, E, Γ)
    E11 = p1.*λ1
    E12 = p2.*λ2
    return sum(m1*E11*Γ[1] + m2*E12*Γ[2])
end

function agg_consumption(c1, c2, m1, m2, λ1, λ2)
    C1 = c1.*λ1
    C2 = c2.*λ2
    return sum(m1*C1 + m2*C2)
end

function agg_savings(s1, s2, m1, m2, λ1, λ2)
    S1 = s1.*λ1
    S2 = s2.*λ2
    return sum(m1*S1 + m2*S2)
end

function simulate_sectoral_shock(r, w1, u, w2lb, w2ub, n_grid, m1, m2, agrid)

    W2 = LinRange(w2lb, w2ub, n_grid)
    WDIST = zeros(n_grid, agrid)
    summary_mat = zeros(12, n_grid)

    # Parameters
    Γ_1 = [1 0.32]
    Γ_2 = [0.45 1]
    β = 0.96
    σ = 2
    σ_ϕ = 1

    Threads.@threads for i in eachindex(W2) # For every wage 2, in parallel

        EV_1, p1_1, p2_1, pu_1, v_1, A, E, C_1, AP_1, λ_1, WD_1 = solver(σ, σ_ϕ, r, w1, W2[i], u, β, Γ_1, agrid); # Solve the Model for type 1
        EV_2, p1_2, p2_2, pu_2, v_2, A, E, C_2, AP_2, λ_2, WD_2 = solver(σ, σ_ϕ, r, w1, W2[i], u, β, Γ_2, agrid); # Solve the model for type 2

        # Compute aggregate measures
        E1 = employment_rate(p1_1, p1_2, m1, m2, λ_1, λ_2, A, E)
        E2 = employment_rate(p2_1, p2_2, m1, m2, λ_1, λ_2, A, E)
        E = E1 + E2
        U =  employment_rate(pu_1, pu_2, m1, m2, λ_1, λ_2, A, E)

        L1 = effective_labor_supply(p1_1, p1_2, m1, m2, λ_1, λ_2, A, E, [Γ_1[1], Γ_2[1]])
        L2 = effective_labor_supply(p2_1, p2_2, m1, m2, λ_1, λ_2, A, E, [Γ_1[2], Γ_2[2]])
        L = L1 + L2
        
        mean_prod_1 = L1/E1
        mean_prod_2 = L2/E2
        prod = L/E

        C = agg_consumption(C_1, C_2, m1, m2, λ_1, λ_2) 
        S = agg_savings(AP_1, AP_2, m1, m2, λ_1, λ_2)



        summary_mat[1, i] = E1
        summary_mat[2, i] = E2
        summary_mat[3, i] = E
        summary_mat[4, i] = U

        summary_mat[5, i] = L1
        summary_mat[6, i] = L2
        summary_mat[7, i] = L

        summary_mat[8, i] = mean_prod_1
        summary_mat[9, i] = mean_prod_2
        summary_mat[10, i] = prod

        summary_mat[11, i] = C
        summary_mat[12, i] = S

        WDIST[i,:] = m1*WD_1 + m2*WD_2
    end

    return summary_mat, W2, WDIST
end


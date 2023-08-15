using Distributed
using JLD2 # extra dependency
addprocs(32; exeflags="--project")

#load survival signatures of subsystems
@load "demo/models/berlin-metro/subsystem-signatures.jld2"

@everywhere begin
    using ResilienceDecisionMaking
    using SurvivalSignature
    using Distributions
    using FiniteDiff

    using LinearAlgebra

    function deleteAdjacency(A_, components)
        A = copy(A_)
        failed = .!components
        A[failed, :] .= Inf
        A[:, failed] .= Inf
        return A
    end

    function efficiency(D)
        n = size(D, 1)

        floyd_warshall!(D)
        D[D .== 0] .= Inf

        D = D .^ -1

        D[diagind(D)] == 0

        return 1 / (n * (n - 1)) * sum(D)
    end

    function floyd_warshall!(D)
        n = size(D, 1)

        @inbounds for k in 1:n, i in 1:n, j in 1:n
            sum = D[i, k] + D[k, j]
            (sum < D[i, j]) && (D[i, j] = sum)
        end

        return D
    end

    ## Overall System
    include("models/berlin-metro/system.jl")

    replace!(A, 0.0 => Inf) # Set non existing edges to infinity for floyd warshall

    ## Subsystem types
    include("models/berlin-metro/types.jl")

    # take reciprocal of adjacency matrix for fastFloyd
    # A = 1 ./ A # !TODO: Where does this come from?

    # Global System Performance Function
    function floydWarshallNetworkEfficiency(A::Matrix{Float64}, comps::BitVector)
        return efficiency(deleteAdjacency(A, comps))
    end

    endowments = [10, 10, 10, 10, 10]

    t = (0.0, 10.0)
    tₙ = 200

    r_max = 22
    α = 0.99

    scale_rate = [0.34, 0.43, 0.36, 0.66]
    Δscale_rate = [0.03, 0.04, 0.034, 0.051]

    function failure_rates(
        times::Vector{T} where {T<:Real}, η::Vector{T} where {T<:Integer}
    )
        param = [
            scale_rate[1] - Δscale_rate[1] * η[1],
            scale_rate[2] - Δscale_rate[2] * η[2],
            scale_rate[3] - Δscale_rate[3] * η[3],
            scale_rate[4] - Δscale_rate[4] * η[4],
        ]
        distributions = Dict(
            1 => Exponential(1 / param[1]),
            2 => Gamma(1 / param[2], 1.2),
            3 => Exponential(1 / param[3]),
            4 => Gamma(1 / param[4], 2.6),
        )

        cumulativeHazard₁(t) = -log(reliability([t], Φ₁, distributions)[1])
        cumulativeHazard₂(t) = -log(reliability([t], Φ₂, distributions)[1])
        cumulativeHazard₃(t) = -log(reliability([t], Φ₃, distributions)[1])
        cumulativeHazard₄(t) = -log(reliability([t], Φ₄, distributions)[1])
        cumulativeHazard₅(t) = -log(reliability([t], Φ₅, distributions)[1])
        cumulativeHazard₆(t) = -log(reliability([t], Φ₆, distributions)[1])

        fdm(f, t) = FiniteDiff.finite_difference_derivative(f, t)

        rates = zeros(size(A, 1), length(times))
        rates[type1, :] .= fdm.(cumulativeHazard₁, times)'
        rates[type2, :] .= fdm.(cumulativeHazard₂, times)'
        rates[type3, :] .= fdm.(cumulativeHazard₃, times)'
        rates[type4, :] .= fdm.(cumulativeHazard₄, times)'
        rates[type5, :] .= fdm.(cumulativeHazard₅, times)'
        rates[type6, :] .= fdm.(cumulativeHazard₆, times)'

        return rates
    end

    function recovery!(components, recovery_times, η)
        for (i, (c, r)) in enumerate(zip(components, recovery_times))
            if ~c
                if r >= r_max - η[5] * 2
                    components[i] = true
                    recovery_times[i] = 1
                else
                    recovery_times[i] = recovery_times[i] + 2
                end
            end
        end

        return nothing
    end
end

res, σ, front = grid_search(
    A,
    endowments,
    t,
    tₙ,
    failure_rates,
    recovery!,
    floydWarshallNetworkEfficiency,
    500,
    α,
    true,
)

function cost_rel(η)
    return 692 * 100 * 1.2^(η[1] - 1) +
           688 * 200 * 1.2^(η[2] - 1) +
           700 * 200 * 1.2^(η[3] - 1) +
           696 * 400 * 1.2^(η[4] - 1)
end

function cost_rec(η)
    return 27 * 100 * 1.2^(η[5] - 1) +
           218 * 200 * 1.2^(η[5] - 1) +
           19 * 300 * 1.2^(η[5] - 1) +
           34 * 400 * 1.2^(η[5] - 1) +
           4 * 500 * 1.2^(η[5] - 1) +
           4 * 600 * 1.2^(η[5] - 1)
end

function cost(η)
    return cost_rel(η) + cost_rec(η)
end

min_cost, min_η = minimal_costs(front, cost)

println("Minimum costs: $min_cost ($min_η)")

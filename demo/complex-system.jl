using ResilienceDecisionMaking

using SurvivalSignature # extra dependency
using Distributions # extra dependency
using FiniteDiff # extra dependency

## Overall System
A = [
    0.0 0.0 1.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
    0.0 0.0 0.0 0.0 1.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
    0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
    0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0
    0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
    0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0
    0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 1.0 0.0 0.0 0.0 0.0
    0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 1.0 0.0
    0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0
    0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0
    0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0
    0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0
    0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0
    0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
]

## Subsystem 1
A1 = zeros(7, 7)
A1[1, [3]] .= 1.0
A1[2, [3]] .= 1.0
A1[3, [4, 5]] .= 1.0
A1[4, [6]] .= 1.0
A1[5, [7]] .= 1.0

types1 = Dict(1 => [2, 3, 4], 2 => [1, 5, 6, 7])

φ = SurvivalSignature.s_t_connectivity([1:7;], [1, 2], [6, 7])
Φ₁ = survivalsignature(A1, types1, φ)

## Subsystem 2
A2 = zeros(7, 7)
A2[1, [3, 4]] .= 1.0
A2[2, [5]] .= 1.0
A2[3, [6]] .= 1.0
A2[4, [6]] .= 1.0
A2[5, [7]] .= 1.0

types2 = Dict(1 => [1, 5, 6], 2 => [2, 3, 4, 7])

φ = SurvivalSignature.s_t_connectivity([1:7;], [1, 2], [6, 7])
Φ₂ = survivalsignature(A2, types2, φ)

## Subsystem 3
A3 = zeros(9, 9)
A3[1, [2, 3]] .= 1.0
A3[2, [4, 5]] .= 1.0
A3[3, [6]] .= 1.0
A3[4, [7]] .= 1.0
A3[5, [8]] .= 1.0
A3[6, [8]] .= 1.0
A3[7, [9]] .= 1.0
A3[8, [9]] .= 1.0

types3 = Dict(1 => [3, 4, 5, 6, 9], 2 => [1, 2, 7, 8])

φ = SurvivalSignature.s_t_connectivity([1:9;], [1], [9])
Φ₃ = survivalsignature(A3, types3, φ)

## Subsystem 4
A4 = zeros(8, 8)
A4[1, [2, 3, 4]] .= 1.0
A4[2, [5, 6]] .= 1.0
A4[3, [5, 6]] .= 1.0
A4[4, [7]] .= 1.0
A4[5, [8]] .= 1.0
A4[6, [8]] .= 1.0
A4[7, [8]] .= 1.0

types4 = Dict(1 => [1, 4, 5, 7], 2 => [2, 3, 6, 8])

φ = SurvivalSignature.s_t_connectivity([1:8;], [1], [8])
Φ₄ = survivalsignature(A4, types4, φ)

## Subsystem 5
A5 = zeros(10, 10)
A5[1, [2]] .= 1.0
A5[2, [3, 4, 5, 6]] .= 1.0
A5[3, [7, 8, 9]] .= 1.0
A5[4, [7, 8, 9]] .= 1.0
A5[5, [7, 8, 9]] .= 1.0
A5[6, [7, 8, 9]] .= 1.0
A5[9, [10]] .= 1.0

types5 = Dict(1 => [2, 3, 4, 5, 10], 2 => [1, 6, 7, 8, 9])

φ = SurvivalSignature.s_t_connectivity([1:10;], [1], [7, 8, 10])
Φ₅ = survivalsignature(A5, types5, φ)

## Subsystem 6

A6 = zeros(7, 7)
A6[1, [2, 3, 4]] .= 1.0
A6[2, [5, 6]] .= 1.0
A6[3, [5, 6]] .= 1.0
A6[4, [5, 6]] .= 1.0
A6[5, [7]] .= 1.0
A6[6, [7]] .= 1.0

types6 = Dict(1 => [1, 2, 3, 7], 2 => [4, 5, 6])

φ = SurvivalSignature.s_t_connectivity([1:7;], [1], [7])
Φ₆ = survivalsignature(A6, types6, φ)

endowments = [10, 10]

t = (0.0, 10.0)
tₙ = 200

r_max = 20
α = 0.90

λ = [0.15, 0.24] # Minor mistake in the paper: [0.15, 0.20]
Δλ = [0.014, 0.022] # Minor mistake in the paper: [0.014, 0.019]

function failure_rates(times::Vector{T} where {T<:Real}, η::Vector{T} where {T<:Integer})
    mean = [1 / (λ[1] - Δλ[1] * η[1]), 1 / (λ[2] - Δλ[2] * η[2])]

    distributions = Dict(1 => Exponential(mean[1]), 2 => Exponential(mean[2]))

    cumulativeHazard₁(t) = -log(reliability([t], Φ₁, distributions)[1])
    cumulativeHazard₂(t) = -log(reliability([t], Φ₂, distributions)[1])
    cumulativeHazard₃(t) = -log(reliability([t], Φ₃, distributions)[1])
    cumulativeHazard₄(t) = -log(reliability([t], Φ₄, distributions)[1])
    cumulativeHazard₅(t) = -log(reliability([t], Φ₅, distributions)[1])
    cumulativeHazard₆(t) = -log(reliability([t], Φ₆, distributions)[1])

    fdm(f, t) = FiniteDiff.finite_difference_derivative(f, t)

    rates = zeros(size(A, 1), length(times))
    rates[[1:5;], :] = ones(5) * fdm.(cumulativeHazard₁, times)'
    rates[[6:9;], :] = ones(4) * fdm.(cumulativeHazard₂, times)'
    rates[[10:11;], :] = ones(2) * fdm.(cumulativeHazard₃, times)'
    rates[12, :] = fdm.(cumulativeHazard₄, times)'
    rates[13, :] = fdm.(cumulativeHazard₅, times)'
    rates[14, :] = fdm.(cumulativeHazard₆, times)'

    return rates
end

function recovery!(components, recovery_times, η)
    for (i, (c, r)) in enumerate(zip(components, recovery_times))
        if ~c
            if r >= r_max
                components[i] = true
                recovery_times[i] = 1
            else
                recovery_times[i] = recovery_times[i] + 2
            end
        end
    end

    return nothing
end

## Execution

res, σ, front = grid_search(
    A,
    endowments,
    t,
    tₙ,
    failure_rates,
    recovery!,
    ResilienceDecisionMaking.s_t_connectivity([1, 2], [14]),
    500,
    α,
    true,
)

function cost(η)
    return 50 * (1000 * 1.2^(η[1] - 1)) + 56 * (2000 * 1.2^(η[2] - 1))
end

min_cost, min_η = minimal_costs(front, cost)

println("Minimum costs: $min_cost ($min_η)")

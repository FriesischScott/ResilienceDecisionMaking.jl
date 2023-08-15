using ResilienceDecisionMaking

A = [
    0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0
    0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0
    0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0
    0.0 0.0 0.0 0.0 1.0 1.0 1.0 1.0
    0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
    0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
    0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
    0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
]

endowments = [20, 20, 20, 20]

t = (0.0, 10.0)
tₙ = 200

λ = 0.8

function failure_rates(times::Vector{T} where {T<:Real}, η::Vector{T} where {T<:Integer})
    rates = [
        (λ - 0.03 * η[1]) * ones(2)...,
        (λ - 0.03 * η[2]),
        (λ - 0.03 * η[3]),
        (λ - 0.03 * η[4]) * ones(4)...,
    ]
    return hcat([rates for _ in 1:length(times)]...)
end

r_max = 21

function recovery!(components, recovery_times, η)
    for (i, (c, r)) in enumerate(zip(components, recovery_times))
        if ~c
            if r >= r_max - 11
                components[i] = true
                recovery_times[i] = 1
            else
                recovery_times[i] = recovery_times[i] + 1
            end
        end
    end

    return nothing
end

α = 0.85

res, σ, front = grid_search(
    A,
    endowments,
    t,
    tₙ,
    failure_rates,
    recovery!,
    ResilienceDecisionMaking.s_t_connectivity([1, 3], [5, 6, 7, 8]),
    500,
    α,
    true,
)

function cost_rel(η)
    return 2 * (800 * 1.2^(η[1] - 1)) +
           800 * 1.2^(η[2] - 1) +
           800 * 1.2^(η[3] - 1) +
           4 * (500 * 1.2^(η[4] - 1))
end

function cost_rec(η)
    return 8 * 600 * 1.2^(11 - 1)
end

function cost(η)
    return cost_rel(η) + cost_rec(η)
end

min_η, min_cost = minimal_costs(front, cost)

println("Minimum costs: $min_cost ($min_η)")

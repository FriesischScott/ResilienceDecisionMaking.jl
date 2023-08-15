# ResilienceDecisionMaking.jl

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7034998.svg)](https://doi.org/10.5281/zenodo.7034998)

This package contains the code and examples for the resilience decision-making method presented in [Salomon et al. (2022)](https://doi.org/10.1016/j.rcns.2022.10.005).

## Installation

The package is registered in the Julia General registry and can be installed from Pkg with

```julia
] add ResilienceDecisionMaking
```

## Getting started

The main interface of the resilience decision making method is the `grid_search` method which returns the computed resilience values as a `Dict` of the endowment configurations, the coefficient of variation and the pareto front. To explain we use the example of the multistage high-speed axial compressor as presented in the paper. First, we define a system representation in this case the adjacency matrix of the system as

```julia
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
```

We further define the maximum values of the endowment configurations we want to explore and the time as a tuple of minimum and maximum as well as the number of time steps

```julia
endowments = [20, 20, 20, 20]

t = (0.0, 10.0)
tₙ = 200
```

Next, we define a function `failure_rates` which returns the current failure rates of all components for all time steps `times` for the current endowment configuration `η`

```julia
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
```

Finally, we define the recovery behaviour of the components through the `recovery!` function. This function must update the `components` and `recovery_times` vectors in place (!) and return `nothing`. Note, that in this example recovery is independent of `η`. See the other `demo` example for how to use `η` to influence the recovery of components.

```julia
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
```

Finally, we can run the `grid_search` where configurations are accepted if the resilience value is above `α = 0.85`. We also pass a structure function evaluating of the system is currently working (s-t-connectivity) and the number of samples (500) to use per resilience computation.

```julia
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
    true, # only return the "dominant" entries of the pareto front
)
```

## References

Salomon, J., Behrensdorf, J., Winnewisser, N., Broggi, M., & Beer, M. (2022). Multidimensional resilience decision-making for complex and substructured systems. Resilient Cities and Structures, 1(3), 61-78, [10.1016/j.rcns.2022.10.005](https://doi.org/10.1016/j.rcns.2022.10.005).

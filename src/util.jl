function minimal_costs(front::Vector{Vector{Int64}}, cost_function::Function)
    min_cost, idx = findmin(cost_function, front)

    return front[idx], min_cost
end

function grid_search(
    A::Matrix,
    η_max::Vector,
    t::Tuple{Real,Real},
    tₙ::Real,
    λ::Function,
    recovery::Function,
    φ::Function,
    samples::Integer,
    α::Float64,
    dominant::Bool=false,
)
    res = Dict{Vector{Int64},Float64}()
    std = Dict{Vector{Int64},Float64}()

    η = ones(Int64, size(η_max))

    progress = ProgressUnknown("Grid search...")
    while true
        r, σ = resilience(A, η, t, tₙ, λ, recovery, φ, samples)

        res[η] = r
        std[η] = σ

        next!(progress)
        if r >= α
            if all(η .== 1)
                return res, std, [η]
            end

            break
        end

        η = η .+ 1

        # No accepted combination can be found
        if any(η .> η_max)
            return res, std, []
        end
    end

    offsets = _generate_offsets(length(η))

    accepted = [η]
    rejected = [η .- 1]

    candidates = [map(o -> (η .- 1) .+ o, offsets)..., map(o -> η .- o, offsets)...]

    it = 0
    if dominant
        length_accepted = 1
        length_rejected = 1
    end

    while ~isempty(candidates)
        it += 1
        η = popfirst!(candidates)

        if any(η .< 1) || any(η .> η_max)
            continue
        end

        if haskey(res, η)
            continue
        end

        acc = _isaccepted(η, accepted)
        rej = false
        if !acc
            rej = _isrejected(η, rejected)
        end

        if !(acc || rej)
            r, σ = resilience(A, η, t, tₙ, λ, recovery, φ, samples)

            res[η] = r
            std[η] = σ

            if r >= α
                push!(accepted, η)
                if !_isrejected(η .- 1, rejected)
                    push!(rejected, η .- 1)
                end
                _add_candidates(candidates, η .- 1, η, η_max, offsets)
            else
                push!(rejected, η)
                if !_isaccepted(η .+ 1, accepted)
                    push!(accepted, η .+ 1)
                end
                _add_candidates(candidates, η, η .+ 1, η_max, offsets)
            end

            next!(progress)
        else
            if acc && (η ∉ accepted && η ∉ rejected)
                if !_isrejected(η .- 1, rejected)
                    push!(rejected, η .- 1)
                end
                if !dominant
                    push!(accepted, η)
                    _add_candidates(candidates, η .- 1, η, η_max, offsets)
                end
            elseif (η ∉ accepted && η ∉ rejected)
                if !_isaccepted(η .+ 1, accepted)
                    push!(accepted, η .+ 1)
                end
                if !dominant
                    push!(rejected, η)
                    _add_candidates(candidates, η, η .+ 1, η_max, offsets)
                end
            end
        end

        if dominant && mod(it, 1000) == 0 # TODO: Make adaptive based on search space
            accepted = [
                accepted[1:length_accepted]...,
                filter(
                    x -> findfirst(y -> y !== x && all(y .<= x), accepted) === nothing,
                    accepted[(length_accepted + 1):end],
                )...,
            ]
            rejected = [
                rejected[1:length_rejected]...,
                filter(
                    x -> findfirst(y -> y !== x && all(y .>= x), rejected) === nothing,
                    rejected[(length_rejected + 1):end],
                )...,
            ]
            length_accepted = length(accepted)
            length_rejected = length(rejected)
        end
    end

    finish!(progress)

    front = [accepted..., map(x -> x .+ 1, rejected)...]
    front = filter(x -> !any(x .> η_max), front)
    if dominant
        front = filter(
            x -> findfirst(y -> y !== x && all(y .<= x), front) === nothing, front
        )
    else
        front = [Set(front)...]
    end

    return res, std, front
end

function _isaccepted(η, accepted)
    return findfirst(x -> all(η .>= x), accepted) !== nothing
end

function _isrejected(η, rejected)
    return findfirst(x -> all(η .<= x), rejected) !== nothing
end

function _add_candidates(candidates, η_low, η_high, η_max, offsets)
    for o in offsets
        c_low = η_low .+ o
        c_high = η_high .- o

        if all(c_low .>= 1) && all(c_low .<= η_max)
            push!(candidates, c_low)
        end

        if all(c_high .>= 1) && all(c_high .<= η_max)
            push!(candidates, c_high)
        end
    end
end

function _generate_offsets(n::Int)
    return collect(eachrow(Matrix{Int64}(I, n, n)))
end

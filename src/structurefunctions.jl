"""
    s_t_connectivity(source::Vector{Int}, target::Vector{Int})

Returns a function (system::Matrix{Float64}, components::BitVector) which will evaluate to
true if a path from any of the *source* nodes to any of the *target* nodes exists in *system*
with only the *components* in working state.
"""
function s_t_connectivity(source::Vector{Int}, target::Vector{Int})
    return function (system::Matrix{Float64}, components::BitVector)
        A = copy(system)
        s = copy(source)
        t = copy(target)

        @inbounds A[.~components, :] .= 0
        @inbounds A[:, .~components] .= 0

        while true
            if any([x in t for x in s])
                return true
            end

            @inbounds A[:, s] .= 0
            s = findall(vec(sum(A[s, :]; dims=1)) .> 0)

            if length(s) == 0
                return false
            end
        end
    end
end

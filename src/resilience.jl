function resilience(
    A::Matrix,
    η::Vector,
    t::Tuple{Real,Real},
    tₙ::Real,
    λ::Function,
    recovery::Function,
    φ::Function,
    samples::Integer,
)
    Q = zeros(tₙ, samples)

    number_of_components = size(A, 1)
    Δt = (t[2] - t[1]) / (tₙ - 1)
    λₜ = λ([t[1]:Δt:t[2];], η) .* Δt

    Q = pmap(1:samples) do _
        res = zeros(tₙ)
        components = trues(number_of_components)
        updated_components = trues(number_of_components)
        recovery_time = ones(number_of_components)

        compute = true

        failures = rand(number_of_components, tₙ)
        time = ones(Int64, number_of_components)

        for j in 1:tₙ
            if compute
                @inbounds res[j] = φ(A, components)
            else
                @inbounds res[j] = res[j - 1]
            end

            @inbounds updated_components[failures[:, j] .< [
                λₜ[c, time[c]] for c in 1:number_of_components
            ]] .= false
            recovery(updated_components, recovery_time, η) # updates components and recovery_times

            changes = components .!= updated_components
            compute = any(changes)

            # Reset times of recovered components``
            recovered_components = .!components .& updated_components
            time[recovered_components] .= 1

            components[changes] = updated_components[changes]

            time .+= 1
        end

        return mean(res ./ res[1])
    end

    res = mean(Q)
    σ = std(Q)

    return res, σ
end

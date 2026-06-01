"""
Defines the main interface to set the [`TimeStep`](@ref) of the simulator.
"""
module TimeSteps
import ..Constants: gyromagnetic_ratio
import ..Geometries: Internal

"""
    Simulation(timestep=(tortuosity=3e-2, gradient=1e-4, permeability=0.5, surface_relaxation=0.01, transfer_rate=0.01, dwell_time=0.1))

Creates an object controlling the timestep of the MCMR simulation.

It can be set by supplying a named tuple to the `timestep` keyword when creating a `Simulation`.

At any time the timestep is guaranteed to be shorter than:
1. `tortuosity` * `size_scale(geometry)`^2 / D, where `size_scale` is the average size of the obstructions and `D` is the `diffusivity`.
2. timestep greater than `permeability` times 1 / (maximum permeability parameter)^2
3. timestep that would allow surface relaxation rate at single collision to be greater than `surface_relaxation`.
4. timestep that would allow magnetisation transfer rate at single collision to be greater than `transfer_rate`.
5. the minimum dwell time of the bound pool times `dwell_time`.
6. (`gradient` /( D * \\gamma^2 * G^2))^(1//3), where \\gamma is the [`gyromagnetic_ratio`](@ref) and `G` is the current `gradient_strength`.
"""
mutable struct TimeStep
    max_timestep :: Float64
    scaled_gradient_precision :: Float64
end


function TimeStep(; 
    diffusivity, geometry, size_scale=nothing, 
    tortuosity=3e-2, gradient=1e-4, verbose=true,
    permeability=0.5, surface_relaxation=0.01,
    transfer=0.01, dwell_time=0.1
    )
    if iszero(diffusivity)
        return TimeStep(Inf, Inf)
    end
    use_size_scale = isnothing(size_scale) ? Internal.size_scale(geometry) : size_scale
    options = (
            tortuosity * use_size_scale^2 / diffusivity,
            Internal.max_timestep_sticking(geometry, diffusivity, transfer),
            max_timestep_permeability(geometry, permeability),
            max_timestep_surface_relaxation(geometry, surface_relaxation),
            Internal.min_dwell_time(geometry) * dwell_time,
        )
    idx = argmin(options)
    if verbose
        lines = ["# Timestep determination"]
        if isinf(minimum(options))
            push!(lines, "No maximum timestep set.")
            push!(lines, "To set a maximum timestep set `timestep=<number>` in the `Simulation` constructor.")
        elseif idx == 1
            push!(lines, "Maximum timestep set by tortuosity constraint based on size of geometry to $(options[1]) ms.")
            if isnothing(size_scale)
                push!(lines, "Size scale of smallest object in the simulation was automatically determined to be $(use_size_scale) um.")
                push!(lines, "If this value is too small, you can set the size scale explicitly by passing on `size_scale=<new_value>` to the `Simulator` constructor.")
            else
                push!(lines, "Size scale was set to $(use_size_scale) um by the user.")
            end
        elseif idx == 2
            push!(lines, "Maximum timestep set by requirement to get accurate rate of transitions from free to bound spins to $(options[2]) ms.")
            push!(lines, "You can alter the sensitivity to magnetisation transfer by changing the value of `timestep=(transfer=...)` from its current value of $(transfer).")
        elseif idx == 3
            push!(lines, "Maximum timestep set by requirement to get accurate rate of spins through the permeable surfaces to $(options[3]) ms.")
            push!(lines, "You can alter the sensitivity to permeability by changing the value of `timestep=(permeability=...)` from its current value of $(permeability).")
        elseif idx == 4
            push!(lines, "Maximum timestep set by requirement to get accurate transverse signal loss due to surface relaxation to $(options[4]) ms.")
            push!(lines, "You can alter the sensitivity to surface relaxation by changing the value of `timestep=(surface_relaxation=...)` from its current value of $(surface_relaxation).")
        elseif idx == 5
            push!(lines, "Maximum timestep set by requirement to limit the transition from bound to free state of spins to $(options[5]) ms.")
            push!(lines, "You can alter the sensitivity to the bound state dwell time by changing the value of `timestep=(dwell_time=...)` from its current value of $(dwell_time).")
        end
        push!(lines, "The actual timestep will be reduced based on the MR sequence(s).")
        @info join(lines, '\n')
    end
    return TimeStep(minimum(options), gradient / diffusivity)
end

(ts::TimeStep)(gradient_strength) = min(
    ts.max_timestep,
    (ts.scaled_gradient_precision / gradient_strength^2)^(1//3)
)

"""
    max_timestep_permeability(geometry, scaling)

Compute the maximum timestep necessary to keep the probability of a spin passing through sufficiently low.
"""
max_timestep_permeability(geometry, scaling) = scaling * Internal.max_permeability_non_inf(geometry)^(-2)

"""
    max_timestep_surface_relaxation(geometry, scaling)

Compute the maximum timestep necessary to maintain accuracy in estimating the transverse signal loss due to surface relaxation.
"""
max_timestep_surface_relaxation(geometry, scaling) = scaling * Internal.max_surface_relaxation(geometry)^(-2)

end
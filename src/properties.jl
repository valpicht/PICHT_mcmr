"""
Basic interface for setting the simulation properties

For the main interface see [`Simulation`](@ref MCMRSimulator.Simulations.Simulation).

Types:
- [`MCMRSimulator.GlobalProperties`](@ref MCMRSimulator.Properties.GlobalProperties)

Methods:
- `R1`
- `R2`
- `off_resonance`
- `stick_probability`
"""
module Properties

"""
    GlobalProperties(; R1=0, R2=0, off_resonance=0)

Stores global MRI and collision properties.
"""
struct GlobalProperties
    R1 :: Float64
    R2 :: Float64
    off_resonance :: Float64
end

function GlobalProperties(; R1=0, R2=0, off_resonance=0) 
    GlobalProperties(R1, R2, off_resonance)
end

for symbol in propertynames(GlobalProperties())
    @eval begin
        $symbol(global_properties::GlobalProperties) = getproperty(global_properties, $(QuoteNode(symbol)))
    end
end

function Base.show(io::IO, prop::GlobalProperties)
    print(io, "GlobalProperties(")
    for (func, units) in [
        (R1, "kHz"),
        (R2, "kHz"),
        (off_resonance, "kHz"),
    ]
        value = func(prop)
        if !iszero(value) && !isinf(value) && !isnan(value)
            print(io, "$(nameof(func))=$(value)$(units), ")
        end
    end
    print(io, ")")
end

"""
    stick_probability(surface_density, dwell_time, diffusivity, timestep)

Computes the probability of a spin getting stuck at the surface given
a [`surface_density`](@ref) and [`dwell_time`](@ref)
as well as the diffusivity (in um^2/ms) and the timestep (in ms).
"""
function stick_probability(surface_density::Number, dwell_time::Number, diffusivity::Number, timestep::Number)
    if iszero(surface_density)
        return 0.
    elseif iszero(dwell_time)
        error("Cannot have a dwell time of zero for any surface containing bound spins.")
    else
        return 1 - exp(-sqrt(π * timestep / diffusivity) * surface_density/ dwell_time / 2)
    end
end
end
"""
Methods to extract properties from the `FixedGeometry`.

MRI properties (using `MRIProperties`):
- `R1`
- `R2`
- `off_resonance`

Collision properties:
- `permeability`
- `surface_relaxation`
- `surface_density`
- `dwell_time`

Other functions:
- `max_timestep_sticking`
"""
module Properties

import StaticArrays: SVector
import ....Properties: GlobalProperties, R1, R2, off_resonance, stick_probability
import ..FixedObstructionGroups: FixedGeometry, FixedObstructionGroup, isinside
import ..Intersections: Intersection
import ..Reflections: Reflection, empty_reflection

"""
    MRIProperties(R1, R2, off_resonance)
    MRIProperties(full_geometry, inside_geometry, global_properties, position, reflection)

Determines the R1, R2, and off-resonance of a particle at given `position` and bound at `Reflection`.
"""
mutable struct MRIProperties
    R1 :: Float64
    R2 :: Float64
    off_resonance :: Float64
end

function MRIProperties(full_geometry::FixedGeometry, inside_geometry::FixedGeometry, glob::GlobalProperties, position::SVector{3, Float64}, stuck_to::Reflection)
    res = MRIProperties(glob.R1, glob.R2, glob.off_resonance)

    # first check stuck properties
    if ~iszero(stuck_to.geometry_index)
        group = full_geometry[stuck_to.geometry_index]
        res.R1 += get_value(Val(:R1), group.surface, stuck_to.obstruction_index)
        res.R2 += get_value(Val(:R2), group.surface, stuck_to.obstruction_index)
        res.off_resonance += get_value(Val(:off_resonance), group.surface, stuck_to.obstruction_index)
    end

    for group in inside_geometry
        props = group.volume
        indices = isinside(group, position, stuck_to)
        for index in indices
            res.R1 += get_value(Val(:R1), props, index)
            res.R2 += get_value(Val(:R2), props, index)
            res.off_resonance += get_value(Val(:off_resonance), props, index)
        end
    end
    return res
end

function get_value(::Val{S}, properties::NamedTuple, index::Integer) where {S}
    res = getproperty(properties, S)
    if res isa Vector
        return res[index]
    else
        return res
    end
end

for symbol in (:R1, :R2, :off_resonance)
    @eval begin
        function $symbol(position::AbstractVector, geometry::FixedGeometry, properties::GlobalProperties=GlobalProperties(), stuck_to::Reflection=empty_reflection)
            MRIProperties(geometry, geometry, properties, SVector{3, Float64}(position), stuck_to).$symbol
        end
    end
end


for symbol in (:permeability, :surface_relaxation, :dwell_time, :surface_density)
    @eval begin
        """
            $($symbol)(geometry, reflection)
        
        Returns the $($symbol) experienced by the spin hitting the surface represented by a `Reflection`.
        The `geometry` has to be a [`FixedGeometry`](@ref).
        """
        function $symbol(group::FixedObstructionGroup, has_hit::Int)
            value = group.surface.$symbol
            if value isa Vector
                return value[has_hit]
            else
                return value
            end
        end

        function $symbol(geometry::FixedGeometry, has_hit::Tuple{Int, Int})
            $symbol(geometry[has_hit[1]], has_hit[2])
        end

        function $symbol(geometry::FixedGeometry, has_hit::Union{Intersection, Reflection})
            $symbol(geometry, (has_hit.geometry_index, has_hit.obstruction_index))
        end
    end
end

"""
    max_timestep_sticking(geometry, diffusivity, scaling)

Returns the maximum timestep that can be used while keeping log(1-[`stick_probability`](@ref)) lower than `scaling`
"""
function max_timestep_sticking(geometry::FixedGeometry, diffusivity::Number, scaling)
    minimum([max_timestep_sticking(group, diffusivity, scaling) for group in geometry])
end

max_timestep_sticking(geometry::FixedGeometry{0}, diffusivity::Number, scaling) = Inf

function max_timestep_sticking(group::FixedObstructionGroup, diffusivity::Number, scaling)
    return scaling * maximum(@. log(1 - stick_probability(group.surface.surface_density, group.surface.dwell_time, diffusivity, 1)))^-2
end

"""
    max_permeability_non_inf(geometry)

Returns the largest permeability within the geometry that is not infinite.
"""
function max_permeability_non_inf(geometry::FixedGeometry)
    maximum(max_permeability_non_inf.(geometry))
end

max_permeability_non_inf(geometry::FixedGeometry{0}) = 0.

function max_permeability_non_inf(group::FixedObstructionGroup)
    if group.surface.permeability isa AbstractVector
        return maximum(filter(a -> ~isinf(a), group.surface.permeability))
    else
        return isinf(group.surface.permeability) ? 0. : group.surface.permeability
    end
end

"""
    max_surface_relaxation(geometry)

Returns the largest surface relaxation parameter within the geometry.
"""
function max_surface_relaxation(geometry::FixedGeometry)
    maximum(max_surface_relaxation.(geometry))
end

max_surface_relaxation(geometry::FixedGeometry{0}) = 0.

function max_surface_relaxation(group::FixedObstructionGroup)
    maximum(group.surface.surface_relaxation)
end


"""
    min_dwell_time(geometry)

Returns the minimum dwell time ignoring any parts of the geometry where there is no bound pool.
"""
function min_dwell_time(geometry::FixedGeometry)
    minimum(min_dwell_time.(geometry))
end

min_dwell_time(geometry::FixedGeometry{0}) = Inf

function min_dwell_time(group::FixedObstructionGroup)
    if group.surface.dwell_time isa Number
        no_bound_pool = (
            group.surface.surface_density isa Number ?
            iszero(group.surface.surface_density) :
            all(iszero.(group.surface.surface_density))
        )
        if no_bound_pool
            return Inf
        else
            return group.surface.dwell_time
        end
    else
        if group.surface.surface_density isa Number
            if iszero(group.surface.surface_density)
                return Inf
            else
                return minimum(group.surface.dwell_time)
            end
        else
            minimum(
                group.surface.dwell_time[group.surface.surface_density .> 0];
                init=Inf
            )
        end
    end
end

end
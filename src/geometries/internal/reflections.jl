"""
Computes the reflection of the spin off some geometry.

This is also used internally to represent a bound spin.

Types:
- `Reflection`

Methods:
- `direction`
- `has_hit`
- `previous_hit`
"""
module Reflections

import StaticArrays: SVector
import LinearAlgebra: ⋅, norm
import ..Intersections: Intersection, has_intersection
import ..Obstructions: detect_intersection

"""
    Reflection(intersection, direction, timestep, ratio_displaced, time_moved, distance_moved, permeable=false)

Represents a reflection of a particle after colliding.

Used internally to represent bound spins.

Parameters:
- intersection: [`Intersection`](@ref) represents the intersection of the trajectory with an obstruction. This will supply:
    - `geometry_index`: index of the geometry that got hit.
    - `obstruction_index`: index of the obstruction that got hit within the geometry.
    - `inside`: whether the inside of the obstruction was hit.
- `direction`: this is stored in an updated form if the obstruction is not permeable (i.e., we have an actual reflection) as indicated by the `permeable` parameter.
- `ratio_displaced`: how large the displacement drawn from the Gaussian distribution was relative to that expected for the timestep.
    This is calculated as "step_size / sqrt(timestep)".
- `time_moved`: how long the particle has already been moving since the last time it was free at the start of a timestep (in ms).
- `distance_moved`: how far the particle has already moved (in um).
"""
struct Reflection
    geometry_index :: Int
    obstruction_index :: Int
    inside :: Bool
    direction :: SVector{3, Float64}
    ratio_displaced :: Float64
    time_moved :: Float64
    distance_moved :: Float64
end

function Reflection(collision::Intersection, direction::SVector{3}, ratio_displaced, time_moved, distance_moved, permeable=false)
    if permeable
        reflection = direction
        inside = ~collision.inside
    else
        reflection = - 2 * (collision.normal ⋅ direction) * collision.normal / norm(collision.normal) ^ 2 .+ direction
        inside = collision.inside
    end
    return Reflection(
        collision.geometry_index, collision.obstruction_index, inside, 
        reflection / norm(reflection), ratio_displaced, time_moved, distance_moved
    )
end

"""
    Reflection(ratio_displaced)

Creates a "virtual" Reflection that represents the movement of a free particle as drawn at the start of a timestep.
"""
Reflection(ratio_displaced) = Reflection(0, 0, false, zero(SVector{3, Float64}), ratio_displaced, 0, 0)

"""
    has_intersection(reflection)

Returns true if this reflection is not empty, i.e., it represents a real collision with an obstruction.
"""
function has_intersection(reflection::Reflection) 
    (i1, i2) = has_hit(reflection)
    return ~(iszero(i1) || iszero(i2))
end

"""
    has_hit(reflection)

Returns a tuple with the indices of the obstruction that the particle is reflecting off.
"""
has_hit(r::Reflection) = (r.geometry_index, r.obstruction_index)

"""
    previous_hit(reflection)

Returns a tuple with the details of the previous reflection relevant for [`detect_intersection`](@ref).
"""
previous_hit(r::Reflection) = (r.geometry_index, r.obstruction_index, r.inside)

"""
    direction(reflection, remaining_time)

Returns the displacement that the spin should end up if it has `remaining_time` ms left in this timestep.
"""
function direction(reflection::Reflection, new_time, diffusivity)
    displacement_size = reflection.ratio_displaced * sqrt(2 * diffusivity * (reflection.time_moved + new_time)) - reflection.distance_moved
    @assert displacement_size > 0
    return reflection.direction * displacement_size
end

const empty_reflection = Reflection(0, 0, false, zero(SVector{3, Float64}), 0., 0., 0.)
end
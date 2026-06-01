"""
All base obstructions are assumed to be centered on the origin (except for triangles in the mesh).

This module defines the `Shift` from this origin to the actual position in space.
"""
module Shifts

import StaticArrays: SVector
import ..FixedObstructions: FixedObstruction, BoundingBox, has_inside, isinside, detect_intersection, radius, size_scale, obstruction_type, random_surface_positions
import ..ObstructionIntersections: ObstructionIntersection

"""
Shifts the obstruction `base` from the origin to `shift`.

The resulting object itself behaves like the base obstruction.
"""
struct Shift{N, O<:FixedObstruction{N}} <: FixedObstruction{N}
    base :: O
    shift :: SVector{N, Float64}
end

Shift(base::FixedObstruction{1}, shift::Float64) = Shift(base, SVector{1}([shift]))
Shift(base::FixedObstruction{N}, shift::AbstractVector) where {N} = Shift(base, SVector{N}(shift))

radius(s::Shift) = radius(s.base)
size_scale(s::Shift) = size_scale(s.base)
function BoundingBox(s::Shift) 
    bb = BoundingBox(s.base)
    return BoundingBox(bb.lower + s.shift, bb.upper + s.shift)
end
has_inside(::Type{Shift{N, O}}) where {N, O} = has_inside(O)
isinside(s::Shift{N}, p::SVector{N}) where {N} = isinside(s.base, p - s.shift)

obstruction_type(::Type{<:Shift{N, O}}) where {N, O} = O

function detect_intersection(s::Shift{N}, start::SVector{N}, dest::SVector{N}) where {N}
    detect_intersection(s.base, start - s.shift, dest - s.shift)
end

function detect_intersection(s::Shift{N}, start::SVector{N}, dest::SVector{N}, inside::Bool) where {N}
    detect_intersection(s.base, start - s.shift, dest - s.shift, inside)
end

function get_quadrant(shift::Shift{N}, repeats::SVector{N}) where {N}
    function f(s, r)
        @assert abs(s) <= r/2
        if s < -r/4
            return -1
        elseif s > r/4
            return 1
        else
            return 0
        end
    end
    SVector{N, Int8}(f.(shift.shift, repeats))
end

function random_surface_positions(s::Shift, density::Number)
    (positions, normals) = random_surface_positions(s.base, density)
    return ([p .+ s.shift for p in positions], normals)
end

end
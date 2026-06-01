"""
Infinite walls and how to intersect with them.
"""
module Walls

import StaticArrays: SVector
import Distributions: Poisson
import Random: rand
import ..FixedObstructions: FixedObstruction, radius, has_inside, isinside, detect_intersection, size_scale, random_surface_positions
import ..ObstructionIntersections: ObstructionIntersection, empty_obstruction_intersections

struct Wall <: FixedObstruction{1}
end

radius(::Wall) = 0.
size_scale(::Wall) = Inf
has_inside(::Type{Wall}) = false
isinside(::Wall, ::SVector{1, Float64}) = false


const negative_normal = SVector{1, Float64}([-1])
const positive_normal = SVector{1, Float64}([1])

function detect_intersection(::Wall, start::AbstractVector, dest::AbstractVector)
    origin = start[1]
    destination = dest[1]

    if origin * destination > 0
        return empty_obstruction_intersections[1]
    end
    inside = origin > 0  # positive side counts as inside for intersections
    ObstructionIntersection(
        abs(origin) / abs(origin - destination),
        inside ? positive_normal : negative_normal,
        inside
    )
end

function detect_intersection(wall::Wall, start::AbstractVector, dest::AbstractVector, ::Bool)
    if abs(start[1]) < 1e-3  
        # if starting within 1 nanometer of this wall, then this is really the same wall, not a repeat
        return empty_obstruction_intersections[1]
    end
    return detect_intersection(wall, start, dest)
end


function random_surface_positions(w::Wall, density::Number)
    nspins = rand(Poisson(density))
    positions = fill(zero(SVector{1, Float64}), nspins)
    normals = fill(positive_normal, nspins)
    return (positions, normals)
end


end
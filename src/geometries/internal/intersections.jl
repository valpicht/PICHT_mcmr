"""
Defines the `Intersection` between a spin path and the geometry.
"""
module Intersections
import StaticArrays: SVector
import ..Obstructions: has_intersection

struct Intersection
    distance :: Float64
    normal :: SVector{3, Float64}
    inside :: Bool
    geometry_index :: Int
    obstruction_index :: Int
end

const empty_intersection = Intersection(Inf, zero(SVector{3, Float64}), false, 0, 0)

has_intersection(intersection::Intersection) = (intersection.distance <= 1) && (intersection.distance >= 0)
end
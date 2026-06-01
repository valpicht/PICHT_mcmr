"""
Defines the intersection between a spin path and an obstruction.

Types:
- `ObstructionIntersection`

Attributes:
- `empty_obstruction_intersections`

Methods:
- `has_intersection`
"""
module ObstructionIntersections
import StaticArrays: SVector
"""
    ObstructionIntersection{N}(distance, normal)

Represents an intersection with a [`FixedObstruction`]{N}.
Distance should be between 0 and 1 indicating where the obstruction intersects.
"""
struct ObstructionIntersection{N}
    distance :: Float64
    normal :: SVector{N, Float64}
    inside :: Bool
end

"""
[`ObstructionIntersection`](@ref) objects that indicate that no intersection was found.

This is a tuple with the ObstructionIntersection for 1, 2, or 3 dimensions.
"""
const empty_obstruction_intersections = (
    ObstructionIntersection(Inf, zero(SVector{1}), false),
    ObstructionIntersection(Inf, zero(SVector{2}), false),
    ObstructionIntersection(Inf, zero(SVector{3}), false)
)

"""
    has_intersection(intersection)

Return true if the intersection object represents a true intersection.
"""
has_intersection(intersection::ObstructionIntersection) = (intersection.distance <= 1) && (intersection.distance >= 0)

end
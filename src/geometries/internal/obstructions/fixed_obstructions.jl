"""
Base types for the API with individually fixed obstructions.

Types:
- `FixedObstruction`

Methods:
- `detect_intersection`
- `has_inside`
- `isinside`
- `radius`
- `BoundingBox`
- `obstruction_type`
- `random_surface_positions`
"""
module FixedObstructions
import StaticArrays: SVector
import ...BoundingBoxes: BoundingBox, isinside

"""
Parent type of all individual obstructions.

N represents the intrinsic dimensionality of the obstruction (e.g., 1 for walls, 2 for cylinders, 3 for spheres/meshes)

Required interface for each obstruction:
- [`detect_intersection`](@ref)(obstruction, start, dest) -> [`ObstructionIntersection`](@ref)
- [`has_inside`](@ref)(obstruction_type) -> Boolean
- [`isinside`](@ref)(obstruction, position) -> Boolean
- At least one of these two:
    - [`radius`](@ref)(obstruction): returns the radius of the obstruc
    - [`BoundingBox`](@ref)(obstruction): returns a Bounding box surrounding the obstruction
"""
abstract type FixedObstruction{N} end
import ..ObstructionIntersections: ObstructionIntersection


"""
    radius(cylinder/sphere...)

Returns the radius of an obstruction.
"""
function radius end

"""
    BoundingBox(obstruction[, vertices])

Returns the smallest bounding box that fully encloses the [`FixedObstruction`](@ref).

`vertices` will only be provided for the [`IndexTriangle`] type.
"""
BoundingBox(o::FixedObstruction{N}) where {N} = BoundingBox{N}(radius(o))
BoundingBox(o::FixedObstruction, ::NamedTuple) = BoundingBox(o)

"""
    detect_intersection(obstruction, start, dest[, inside])

Determines whether the line from `start` to `dest` intersects with `obstruction`.
If no intersection is found, [`empty_obstruction_intersection`](@ref) is returned.
Otherwise a [`ObstructionIntersection`](@ref) object is returned.

If `inside` is set, the last collision in the trajectory was with this `obstruction`.
`inside` indicates whether this previous collision was with the inside of the `obstruction`.
This is passed on to prevent the particle from colliding with the same obstruction multiple times without any intermediate movement.
"""
detect_intersection(obstruction::FixedObstruction{N}, start::SVector{N}, dest::SVector{N}, ::NamedTuple, prev_inside::Union{Nothing, Bool}) where{N} = detect_intersection(obstruction, start, dest, prev_inside)
detect_intersection(obstruction::FixedObstruction{N}, start::SVector{N}, dest::SVector{N}, ::NamedTuple) where{N} = detect_intersection(obstruction, start, dest)

"""
    has_inside(obstruction_type)

Returns true if a particular sub-type of [`FixedObstruction`](@ref) has an inside.
"""
function has_inside end

"""
    isinside(obstruction, position)

Returns true if the `position` is inside the `obstruction`.
This will be inaccurate if the particle with that position is stuck on the surface of the obstruction.
"""
isinside(o::FixedObstruction{N}, pos::AbstractVector, args...) where {N} = isinside(o, SVector{N}(pos), args...)
isinside(o::FixedObstruction{N}, pos::SVector{N}, ::NamedTuple) where {N} = isinside(o, SVector{N}(pos))

"""
    obstruction_type(obstruction/geometry)

Returns the type of the obstruction that this geometry is built up out of
"""
obstruction_type(::Type{O}) where {O <: FixedObstruction} = O

"""
    size_scale(obstruction/geometry)

Returns the size scale of the smallest object within the geometry (typically the radius or the repeating length scale).
"""
size_scale(o::FixedObstruction) = radius(o)


"""
    random_surface_positions(obstruction[, args], surface_density)

Returns a vector of positions and normals on the normals of a specific obstruction.
Normals should always assume that the particle is on the inside.
"""
random_surface_positions(obstruction::FixedObstruction, ::NamedTuple, surface_density::Number) = random_surface_positions(obstruction, surface_density)


end
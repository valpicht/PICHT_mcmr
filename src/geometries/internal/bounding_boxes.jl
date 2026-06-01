"""
Defines the bounding box and its interface.

Types:
- [`BoundingBox`](@ref): N-dimensional box fully containing an obstruction.

Functions
- [`MCMRSimulator.isinside`](@ref MCMRSimulator.Geometries.Internal.BoundingBoxes.isinside): returns whether point/spin is inside a bounding box (or other obstruction)
- `could_intersect`: returns whether spin trajectory could intersect with obstruction in bounding box
"""
module BoundingBoxes
import StaticArrays: SVector

"""
    BoundingBox(lower::Vector, upper::Vector)
    BoundingBox([center::Vector, ]radius::Number)

Represents a box in N-dimensional space.

Check whether particles are inside using [`isinside`](@ref).

This can be used to define where spins should be generated in a [`Snapshot`](@ref MCMRSimulator.Spins.Snapshot) (see [`readout`](@ref MCMRSimulator.Evolve.readout)).
"""
struct BoundingBox{N}
    lower :: SVector{N, Float64}
    upper :: SVector{N, Float64}
    function BoundingBox(lower::SVector{N, Float64}, upper::SVector{N, Float64}) where {N}
        new{N}(lower, upper)
    end
end

function BoundingBox(lower::AbstractVector, upper::AbstractVector)
    ndim = length(lower)
    BoundingBox(SVector{ndim, Float64}(lower), SVector{ndim, Float64}(upper))
end

function BoundingBox(center::AbstractVector, radius::Number)
    BoundingBox(center .- radius, center .+ radius)
end

function BoundingBox{N}(radius::Number) where {N}
    BoundingBox(fill(-radius, N), fill(radius, N))
end
BoundingBox(radius::Number) = BoundingBox{3}(radius)

BoundingBox(bb::BoundingBox) = bb

lower(bb::BoundingBox{N}) where {N} = bb.lower
upper(bb::BoundingBox{N}) where {N} = bb.upper

function BoundingBox(bounding_boxes::Union{AbstractVector{<:BoundingBox{M}}, NTuple{N, BoundingBox{M}}}) where {N, M}
    lowers = map(lower, bounding_boxes)
    uppers = map(upper, bounding_boxes)
    BoundingBox(
        SVector{M, Float64}(minimum([l[i] for l in lowers]) for i in 1:M), 
        SVector{M, Float64}(maximum([u[i] for u in uppers]) for i in 1:M)
    )
end

isinside(bb::BoundingBox, pos::AbstractVector) = all(pos .>= lower(bb)) && all(pos .<= upper(bb))

"""
    could_intersect(bounding_box, start, dest)

Returns true if a line connecting `start` to `dest` could intersect with an obstruction in `bounding_box`.
This function only does some very fast, high-level checks.
Just because it returns true does not mean that there is an interesection of the obstruction (or even the bounding box) with the line.
However, if it returns false, there is guaranteed to be no intersection.
"""
function could_intersect(bb::BoundingBox{1}, start::AbstractVector, dest::AbstractVector)
    return (
        (start[1] <= bb.upper[1] || dest[1] <= bb.upper[1]) &&
        (start[1] >= bb.lower[1] || dest[1] >= bb.lower[1])
    )
end

function could_intersect(bb::BoundingBox{2}, start::AbstractVector, dest::AbstractVector)
    return (
        (start[1] <= bb.upper[1] || dest[1] <= bb.upper[1]) &&
        (start[1] >= bb.lower[1] || dest[1] >= bb.lower[1]) &&
        (start[2] <= bb.upper[2] || dest[2] <= bb.upper[2]) &&
        (start[2] >= bb.lower[2] || dest[2] >= bb.lower[2])
    )
end

function could_intersect(bb::BoundingBox{3}, start::AbstractVector, dest::AbstractVector)
    return (
        (start[1] <= bb.upper[1] || dest[1] <= bb.upper[1]) &&
        (start[1] >= bb.lower[1] || dest[1] >= bb.lower[1]) &&
        (start[2] <= bb.upper[2] || dest[2] <= bb.upper[2]) &&
        (start[2] >= bb.lower[2] || dest[2] >= bb.lower[2]) &&
        (start[3] <= bb.upper[3] || dest[3] <= bb.upper[3]) &&
        (start[3] >= bb.lower[3] || dest[3] >= bb.lower[3])
    )
end

end
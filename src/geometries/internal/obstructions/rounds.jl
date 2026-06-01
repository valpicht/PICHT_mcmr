"""
Defines cylinders and spheres and how to intersect with them.
"""
module Rounds

import StaticArrays: SVector
import Distributions: Poisson
import Random: rand
import ..FixedObstructions: FixedObstruction, radius, has_inside, isinside, detect_intersection, random_surface_positions
import ..ObstructionIntersections: ObstructionIntersection, empty_obstruction_intersections

struct Round{N} <: FixedObstruction{N}
    radius :: Float64
end

const Cylinder = Round{2}
const Sphere = Round{3}

radius(r::Round) = r.radius
has_inside(::Type{<:Round}) = true
isinside(r::Round{N}, pos::SVector{N}) where {N} = sum(pos .* pos) < radius(r)^2

function detect_intersection(round::Round{N}, start::SVector{N}, dest::SVector{N}) where {N}
    inside = sum(start .* start) <= radius(round)^2
    return detect_intersection(round, start, dest, inside)
end

function detect_intersection(round::Round{N}, start::SVector{N}, dest::SVector{N}, inside::Bool) where {N}
    diff = dest - start

    # terms for quadratic equation for where distance squared equals radius squared d^2 = a s^2 + b s + c == radius ^ 2
    a = sum(diff .* diff)
    b = sum(2 .* start .* diff)
    c = sum(start .* start)
    determinant = b * b - 4 * a * (c - radius(round) ^ 2)
    if determinant < 0
        return empty_obstruction_intersections[N]
    end
    sd = sqrt(determinant)
    ai = inv(a)

    solution = (inside ? (-b + sd) : (-b - sd)) * 1//2 * ai
    if solution > 1 || solution <= 0
        return empty_obstruction_intersections[N]
    end
    normal = solution * dest + (1 - solution) * start
    return ObstructionIntersection(
        solution,
        inside ? -normal : normal,
        inside
    )
end

function random_surface_positions(c::Cylinder, density::Number)
    surface = 2π * c.radius
    nspins = rand(Poisson(surface * density))
    theta = rand(nspins) .* 2π
    normals = [SVector{2, Float64}([cos(t), sin(t)]) for t in theta]
    positions = normals .* (-c.radius)
    return (positions, normals)
end

function random_on_sphere()
    z = rand(Float64) * 2 - 1
    r = sqrt(1 - z*z)
    theta = rand(Float64) * Float64(2 * π)
    (s, c) = sincos(theta)
    return SVector{3, Float64}(
        r * s,
        r * c,
        z
    )
end

function random_surface_positions(s::Sphere, density::Number)
    surface = 4π * s.radius * s.radius
    nspins = rand(Poisson(surface * density))
    normals = [random_on_sphere() for _ in 1:nspins]
    positions = normals .* (-s.radius)
    return (positions, normals)
end

end
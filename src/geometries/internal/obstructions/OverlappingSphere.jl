# Internal geometry support for substrates represented as unions of overlapping
# spheres. This is mainly used for SWC-derived morphologies approximated by
# sphere centers and radii, without constructing a full mesh.

import StaticArrays: SVector
import LinearAlgebra: norm
import Distributions: Poisson
import Random: rand
import .FixedObstructions: FixedObstruction, radius, has_inside
import .ObstructionIntersections: ObstructionIntersection, empty_obstruction_intersections

"""
    OverlappingSpheres

Internal fixed-geometry representation of a set of overlapping spheres in 3D.

Each sphere is stored as a tuple:

    (center, radius)

where:
- `center` is an `SVector{3, Float64}` containing the sphere center coordinates
- `radius` is the sphere radius

This object is used to represent substrates reconstructed as unions of spheres,
for example when approximating SWC morphologies without building a full mesh.
"""
struct OverlappingSpheres <: FixedObstruction{3}
    spheres :: Vector{Tuple{SVector{3, Float64}, Float64}}
end

"""
    radius(os::OverlappingSpheres)

Return a global size estimate for the full set of overlapping spheres.

The value is computed as the maximum of:

    norm(center) + radius

over all spheres, which gives an outer bound on the geometry extent relative to
the origin.

# Arguments
- `os`: an `OverlappingSpheres` object

# Returns
- A scalar radius enclosing the full sphere union
"""
function radius(os::OverlappingSpheres)
    maximum(norm(pos) + r for (pos, r) in os.spheres)
end

"""
    has_inside(::Type{<:OverlappingSpheres})

Declare that `OverlappingSpheres` defines a volume with an interior region.

This is needed by the simulator so that the geometry is treated as an object
that can contain spins inside it.
"""
has_inside(::Type{<:OverlappingSpheres}) = true

"""
    random_on_sphere_overlapping()

Sample a random unit vector uniformly on the surface of a sphere.
"""
function random_on_sphere_overlapping()
    z = rand(Float64) * 2 - 1
    r = sqrt(1 - z*z)
    theta = rand(Float64) * Float64(2 * π)
    s, c = sincos(theta)

    return SVector{3, Float64}(
        r * s,
        r * c,
        z
    )
end

"""
    FixedObstructions.random_surface_positions(os::OverlappingSpheres, density::Number)

Generate random positions on the exposed surface of a union of overlapping
spheres.

This follows the same interface as `random_surface_positions(::Sphere, density)`:
it returns `(positions, normals)`.

Important:
- candidate points are sampled on the surfaces of all individual spheres;
- points hidden inside another sphere are rejected;
- therefore only the external boundary of the union contributes.
"""
function FixedObstructions.random_surface_positions(os::OverlappingSpheres, density::Number)

    # Match behaviour of normal spheres: if density is zero, return empty vectors.
    if density <= 0 || isempty(os.spheres)
        return (
            SVector{3, Float64}[],
            SVector{3, Float64}[]
        )
    end

    centers = [center for (center, r) in os.spheres]
    radii = [r for (center, r) in os.spheres]

    # Full surface area of each sphere. We sample from the full surfaces,
    # then reject hidden points. This gives the right expected density on
    # the exposed surface by thinning.
    areas = 4π .* radii .* radii
    total_area = sum(areas)

    n_candidates = rand(Poisson(total_area * density))

    positions = SVector{3, Float64}[]
    normals = SVector{3, Float64}[]

    if n_candidates == 0
        return (positions, normals)
    end

    cumulative_areas = cumsum(areas)
    tol = 1e-12

    for _ in 1:n_candidates
        # Choose sphere proportional to surface area.
        u = rand(Float64) * total_area
        i = searchsortedfirst(cumulative_areas, u)

        center = centers[i]
        r = radii[i]

        # Match the convention used in rounds.jl:
        # positions = normals .* (-radius)
        normal = random_on_sphere_overlapping()
        position = center - r * normal

        # Reject points that are hidden inside another sphere.
        exposed = true

        for j in eachindex(centers)
            j == i && continue

            other_center = centers[j]
            other_r = radii[j]

            if sum((position - other_center) .* (position - other_center)) < other_r^2 - tol
                exposed = false
                break
            end
        end

        if exposed
            push!(positions, position)
            push!(normals, normal)
        end
    end

    return (positions, normals)
end

"""
    FixedObstructions.isinside(os::OverlappingSpheres, pos)

Test whether a 3D position lies inside the union of the overlapping spheres.

A point is considered inside if it lies inside at least one sphere.

# Arguments
- `os`: an `OverlappingSpheres` object
- `pos`: 3D position to test

# Returns
- `true` if `pos` is inside at least one sphere
- `false` otherwise
"""
function FixedObstructions.isinside(os::OverlappingSpheres, pos::SVector{3, Float64})
    any(sum((pos - center) .* (pos - center)) < r^2 for (center, r) in os.spheres)
end

"""
    sphere_intersection(center, r, start, dest, inside)

Compute the intersection between a line segment and a single sphere.

The segment goes from `start` to `dest`. The sphere is defined by its `center`
and radius `r`.

The parameter `inside` indicates whether the starting point is already inside
that sphere:
- if `inside == false`, the function looks for an entering intersection
- if `inside == true`, the function looks for an exiting intersection

# Arguments
- `center`: sphere center
- `r`: sphere radius
- `start`: starting point of the segment
- `dest`: end point of the segment
- `inside`: whether the segment starts inside the sphere

# Returns
- `nothing` if there is no valid intersection on the segment
- otherwise a tuple `(solution, hit, normal, inside)` where:
  - `solution` is the segment parameter in `(0, 1]`
  - `hit` is the intersection point
  - `normal` is the surface normal at the hit point
  - `inside` is the input inside/outside state used for this sphere
"""
function sphere_intersection(center::SVector{3, Float64}, r::Float64,
                             start::SVector{3, Float64}, dest::SVector{3, Float64},
                             inside::Bool)

    # Express the segment in sphere-centered coordinates
    start_local = start - center
    dest_local = dest - center
    diff = dest_local - start_local

    # Solve the quadratic equation for line-sphere intersection
    a = sum(diff .* diff)
    b = sum(2 .* start_local .* diff)
    c = sum(start_local .* start_local)

    determinant = b * b - 4 * a * (c - r^2)
    if determinant < 0
        return nothing
    end

    sd = sqrt(determinant)
    ai = inv(a)

    # Choose the appropriate root depending on whether we are entering or exiting
    solution = (inside ? (-b + sd) : (-b - sd)) * 0.5 * ai
    if solution > 1 || solution <= 1e-12
        return nothing
    end

    # Reconstruct the intersection point in global coordinates
    hit = solution * dest + (1 - solution) * start
    normal = hit - center

    # If the particle starts inside the sphere, the relevant surface normal
    # is reversed because the crossing is an exit event.
    return (solution, hit, inside ? -normal : normal, inside)
end

"""
    FixedObstructions.detect_intersection(os, start, dest)

Detect the first relevant intersection between a trajectory segment and the
union of overlapping spheres.

This method first determines whether `start` is inside the union of spheres,
then delegates to the method that takes the union-state explicitly.

# Arguments
- `os`: an `OverlappingSpheres` object
- `start`: starting point of the segment
- `dest`: end point of the segment

# Returns
- An `ObstructionIntersection` if a valid boundary crossing is found
- An empty intersection otherwise
"""
function FixedObstructions.detect_intersection(os::OverlappingSpheres,
                                               start::SVector{3, Float64},
                                               dest::SVector{3, Float64})
    inside_union = FixedObstructions.isinside(os, start)
    return FixedObstructions.detect_intersection(os, start, dest, inside_union)
end

"""
    FixedObstructions.detect_intersection(os, start, dest, inside_union)

Detect the first valid intersection between a trajectory segment and the union
boundary of overlapping spheres.

This method handles a subtle but important case: when spheres overlap, some
parts of individual sphere surfaces are hidden inside the union and should not
count as true entry/exit boundaries.

The algorithm:
1. tests intersections with all individual spheres,
2. keeps only intersections consistent with the current union state,
3. removes intersections hidden inside another sphere,
4. returns the nearest valid intersection.

# Arguments
- `os`: an `OverlappingSpheres` object
- `start`: starting point of the segment
- `dest`: end point of the segment
- `inside_union`: whether the segment starts inside the sphere union

# Returns
- An `ObstructionIntersection` corresponding to the nearest valid crossing
- An empty intersection if no valid crossing exists
"""
function FixedObstructions.detect_intersection(os::OverlappingSpheres,
                                               start::SVector{3, Float64},
                                               dest::SVector{3, Float64},
                                               inside_union::Bool)

    best = nothing
    tol = 1e-12

    for (i, (center, r)) in enumerate(os.spheres)
        # Check whether the start point is inside this specific sphere
        inside_this = sum((start - center) .* (start - center)) <= r^2 + tol
        candidate = sphere_intersection(center, r, start, dest, inside_this)

        if isnothing(candidate)
            continue
        end

        solution, hit, normal, inside_candidate = candidate

        # Ignore near-zero self-intersections that can create artificial
        # repeated collisions or bounce loops.
        if solution <= tol
            continue
        end

        # Keep only intersections that match the union-level transition:
        # - outside the union -> entering hit
        # - inside the union  -> exiting hit
        if inside_candidate != inside_union
            continue
        end

        # If the hit point lies strictly inside another sphere, then this
        # local sphere surface is hidden inside the union and should not be
        # treated as a true boundary crossing.
        hidden_by_other = false
        for (j, (other_center, other_r)) in enumerate(os.spheres)
            if i == j
                continue
            end
            if sum((hit - other_center) .* (hit - other_center)) < other_r^2 - tol
                hidden_by_other = true
                break
            end
        end

        if hidden_by_other
            continue
        end

        # Keep the earliest valid intersection along the segment
        if isnothing(best) || solution < best[1]
            best = (solution, normal, inside_candidate)
        end
    end

    if isnothing(best)
        return empty_obstruction_intersections[3]
    end

    solution, normal, inside = best
    return ObstructionIntersection(solution, normal, inside)
end
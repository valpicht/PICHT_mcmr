"""
Triangle within a mesh and how to intersect with them.
"""
module Triangles

import StaticArrays: SVector
import LinearAlgebra: norm, cross, ⋅
import Distributions: Poisson
import Random: rand
import ..FixedObstructions: FixedObstruction, BoundingBox, has_inside, isinside, detect_intersection, size_scale, random_surface_positions, radius
import ..ObstructionIntersections: ObstructionIntersection, empty_obstruction_intersections

abstract type Triangle <: FixedObstruction{3} end


struct IndexTriangle <: Triangle
    indices :: SVector{3, Int32}
    component :: Int32
end

struct FullTriangle <: Triangle
    a :: SVector{3, Float64}
    b :: SVector{3, Float64}
    c :: SVector{3, Float64}
end

Base.getindex(t::FullTriangle, i::Int) = (t.a, t.b, t.c)[i]

const Vertices = Vector{SVector{3, Float64}}

function FullTriangle(indices::SVector{3, <:Integer}, vertices::Vertices)
    FullTriangle(
        vertices[indices[1]],
        vertices[indices[2]],
        vertices[indices[3]]
    )
end

function FullTriangle(indices_triangle::IndexTriangle, vertices::Vertices)
    FullTriangle(indices_triangle.indices, vertices)
end

function BoundingBox(t::FullTriangle)
    BoundingBox(
        SVector{3}([
            min(t.a[1], t.b[1], t.c[1]),
            min(t.a[2], t.b[2], t.c[2]),
            min(t.a[3], t.b[3], t.c[3]),
        ]),
        SVector{3}([
            max(t.a[1], t.b[1], t.c[1]),
            max(t.a[2], t.b[2], t.c[2]),
            max(t.a[3], t.b[3], t.c[3]),
        ])
    )
end

BoundingBox(t::IndexTriangle, args::NamedTuple) = BoundingBox(FullTriangle(t, args.vertices))
    
function radius(t::FullTriangle) 
    bb = BoundingBox(t)
    return maximum(bb.upper .- bb.lower) / 2
end

has_inside(::Type{<:Triangle}) = true
isinside(::Triangle, ::SVector{3, Float64}) = false

"""
    triangle_size(full_triangle)
    triangle_size(p1, p2, p3)

Computes the size of a [`FullTriangle`](@ref) formed by three points
"""
triangle_size(ft::FullTriangle) = triangle_size(ft.a, ft.b, ft.c)

function triangle_size(p1, p2, p3)
    e1 = p2 - p1
    e2 = p3 - p1
    sz = cross(e1, e2)
    return norm(sz) / 2
end

"""
    normal(p1, p2, p3)

Computes the normal of a triangle formed by the three points
"""
function normal(p1 :: AbstractVector, p2 :: AbstractVector, p3 :: AbstractVector)
    e1 = p2 - p1
    e2 = p3 - p1
    tonorm = cross(e1, e2)
    return tonorm ./ norm(tonorm)
end

function normal(t::FullTriangle)
    return normal(t.a, t.b, t.c)
end


function detect_intersection(triangle::FullTriangle, start::SVector{N}, dest::SVector{N}, inside::Union{Nothing, Bool}) where {N}
    return detect_intersection_partial(triangle, start, dest, inside)[1]
end

function detect_intersection(triangle::IndexTriangle, start::SVector{N}, dest::SVector{N}, args::NamedTuple, inside::Union{Nothing, Bool}=nothing) where {N}
    return detect_intersection(FullTriangle(triangle, args.vertices), start, dest, inside)
end

"""
Computes the intersection for a triangle in a mesh (returned by [`detect_intersection`](@ref)).
This function also returns an additional bool.
This bool will be true if the intersection is exactly at the edge of the triangle.
"""
function detect_intersection_partial(triangle::FullTriangle, start::SVector{N}, dest::SVector{N}, inside=nothing) where {N}
    t_normal = normal(triangle)
    dist_plane = t_normal ⋅ triangle.a

    dist_start = t_normal ⋅ start
    dist_dest = t_normal ⋅ dest
    if abs(dist_start - dist_dest) < 1e-8
        # moving parallel to the triangle
        return (empty_obstruction_intersections[3], false)
    end
    if ~isnothing(inside) && (dist_plane - dist_start) < 1e-3
        # if starting point is within 1 nm of the previous intersection than this is the same triangle, not a repeat
        return (empty_obstruction_intersections[3], false)
    end

    time = (dist_plane - dist_start ) / (dist_dest - dist_start)

    if time < 0 || time > 1
        # does not reach the plane of the triangle
        return (empty_obstruction_intersections[3], false)
    end

    # where the trajectory intersects the plane
    intersect_point = @. time * dest + (1 - time) * start

    partial = false
    for (dim, d1) in (
            (1, 2),
            (2, 3),
            (3, 1),
        )
        edge = triangle[d1] - triangle[dim]
        to_point = intersect_point - triangle[dim]
        along_normal = cross(edge, to_point)
        inpr = along_normal ⋅ t_normal
        if (inpr) < 0
            # intersect point is on the wrong side of this edge and hence not in the triangle
            return (empty_obstruction_intersections[3], false)
        elseif iszero(inpr)
            partial = true
        end
    end
    inside = dist_dest > dist_start
    return (ObstructionIntersection(
        time,
        inside ? -t_normal : t_normal,
        inside
    ), partial)
end

function random_surface_positions(triangle::IndexTriangle, args::NamedTuple, density::Number)
    random_surface_positions(FullTriangle(triangle, args.vertices), density)
end

function random_surface_positions(ft::FullTriangle, density::Number)
    surface = triangle_size(ft.a, ft.b, ft.c)
    nspins = rand(Poisson(surface * density))

    e1 = ft.b .- ft.a
    e2 = ft.c .- ft.a
    function draw_position()
        (u1, u2) = rand(2)
        if u1 + u2 > 1
            u1 = 1 - u1
            u2 = 1 - u2
        end
        return SVector{3}(@. ft.a + u1 * e1 + u2 * e2)
    end
    positions = [draw_position() for _ in 1:nspins]
    normals = fill(-normal(ft.a, ft.b, ft.c), nspins)
    return (positions, normals)
end

"""
    curvature(triangles, vertices[index1, index2])

Computes the curvature between two neighbouring triangles (`index1` and `index2`).
If no indices are provided computes the mean curvature over the whole surface.
"""
function curvature(triangles, vertices, index1, index2)
    pos1 = map(i->vertices[i], triangles[index1])
    pos2 = map(i->vertices[i], triangles[index2])
    pos_offset = @. (pos1[1] + pos1[2] + pos1[3] - pos2[1] - pos2[2] - pos2[3]) / 3
    normal_offset = normal(pos1...) - normal(pos2...)
    return (pos_offset ⋅ normal_offset) / norm(pos_offset)^2
end

function curvature(triangles, vertices) 
    sizes = [triangle_size(vertices[t[1]], vertices[t[2]], vertices[t[3]]) ./ 3 for t in triangles]
    all_curv = [(sizes[i1] + sizes[i2], curvature(triangles, vertices, i1, i2)) for (i1, i2) in neighbours(triangles)]
    return sum([w .* c for (w, c) in all_curv]) / sum([w for (w, _) in all_curv])
end

curvature(triangles::AbstractVector{IndexTriangle}, vertices) = curvature([t.indices for t in triangles], vertices)

"""
    neighbours(triangles)

Return pairs of indices of triangles that share an edge.
"""
function neighbours(triangles::AbstractVector)
    norm_edge((i1, i2)) = i1 > i2 ? (i2, i1) : (i1, i2)
    edges(t) = norm_edge.([(t[1], t[2]), (t[2], t[3]), (t[3], t[1])])

    edges_to_triangle = Dict{Tuple{Int, Int}, Int}()
    pairs = Tuple{Int, Int}[]
    for (i, t) in enumerate(triangles)
        for e in edges(t)
            if e in keys(edges_to_triangle)
                push!(pairs, (edges_to_triangle[e], i))
            else
                edges_to_triangle[e] = i
            end
        end
    end
    return pairs
end


end
"""
Methods that extend the [`Mesh`](@ref) constructor to convert different types into meshes.
"""
module ToMesh
import StaticArrays: SVector, MVector
import BSplineKit: BSplineOrder, interpolate, Derivative
import LinearAlgebra: cross, norm, ⋅
import Statistics: mean
import ....Methods: get_rotation
import ..Obstructions: Mesh, value_as_vector, BendyCylinder, Cylinder, Cylinders, Wall, Walls, Sphere, Spheres

"""
    Mesh(other_obstruction; kwargs...)

Approximates any other obstruction type by a mesh.
Sometimes a sequence of [`Mesh`](@ref) objects will be returned instead
(e.g., one mesh for each sphere/cylinder in [`Spheres`](@ref)/[`Cylinders`](@ref) objects).

## Keyword arguments
For [`Cylinders`](@ref):
- `nsamples`: number of mesh vertices along the circumference (default: 100).
- `height`: height of mesh triangles in um along the long axis (default: average cylinder circumference divided by `nsamples`).

For [`Walls`](@ref):
- `height`: heigh of mesh triangle in um

For [`Spheres`](@ref):
- `nasamples`: resulting mesh will contain at least this many triangles (default: 180).
"""
function Mesh(bendy_cylinder::BendyCylinder)
    control_points = value_as_vector(bendy_cylinder.control_point)
    radii = value_as_vector(bendy_cylinder.radius)
    triangle_size = 2π * mean(radii) / bendy_cylinder.nsamples.value

    nadd = isnothing(bendy_cylinder.closed.value) ? 0 : div(bendy_cylinder.spline_order.value, length(bendy_cylinder.control_point), RoundUp)
    if iszero(nadd)
        full_control_points = SVector{3}.(control_points)
        full_radii = Float64.(radii)
    else
        @assert !isnothing(bendy_cylinder.repeats.value)
        full_control_points = [SVector{3}(cp .+ bendy_cylinder.repeats.value .* bendy_cylinder.closed.value .* index_add) for index_add in -nadd:nadd for cp in control_points]
        full_radii = [r for _ in -nadd:nadd for r in radii]
    end

    dist = [0, cumsum(norm.(full_control_points[2:end] .- full_control_points[1:end-1]))...]
    dist .-= dist[nadd * length(control_points) + 1]  # set first of original control points to reference point with distance of 0

    fpos = interpolate(dist, full_control_points, BSplineOrder(bendy_cylinder.spline_order.value))
    frad = interpolate(dist, full_radii, BSplineOrder(bendy_cylinder.spline_order.value))

    nsamples = bendy_cylinder.nsamples.value
    theta = range(start=0, stop=2π, length=nsamples + 1)[1:end-1]
    if iszero(nadd)
        dist_control = dist
    else
        dist_control = dist[nadd * length(control_points) + 1: (nadd + 1) * length(control_points) + 1]
    end
    dist_eval = [0., [
        d
        for (d1, d2) in zip(dist_control[1:end-1], dist_control[2:end])
        for d in range(start=d1, stop=d2, length=max(Int(div(d2 - d1, triangle_size, RoundUp)), 2))[2:end]
    ]...]

    vertices = SVector{3, Float64}[]
    triangles = SVector{3, Int}[]

    fdirection = Derivative(1) * fpos

    # start off vec_theta0 in arbitrary direction
    vec_theta0_try = cross(fdirection(0), SVector{3}((1., 0., 0.)))    
    if all(iszero.(vec_theta0_try))
        vec_theta0_try = cross(fdirection(0), SVector{3}((0., 1., 0.)))    
    end
    vec_theta0 = MVector{3}(vec_theta0_try ./ norm(vec_theta0_try))

    function get_vertices(dist)
        # adjust vec_theta0 to be orthogonal to current direction
        direction_unnorm = fdirection(dist)
        direction = direction_unnorm ./ norm(direction_unnorm)
        vec_theta0 .-= direction .* (direction ⋅ vec_theta0)
        vec_theta0 ./= norm(vec_theta0)
        centroid = fpos(dist)
        radius_uncorrected = frad(dist)
        # correct radius to preserve final volume
        half_theta_step = π / nsamples
        radius = radius_uncorrected * (half_theta_step / (sin(half_theta_step) * cos(half_theta_step)))
        last_vec = cross(direction, vec_theta0)
        return [SVector{3}((radius * sin(t)) .* last_vec .+ (radius * cos(t)) .* vec_theta0 .+ centroid) for t in theta]
    end

    vertices = vcat(get_vertices.(dist_eval)...)

    single_triangle_set = [
        SVector{3}(nsamples, 1, nsamples * 2),
        [SVector{3}((add, 1 + add, nsamples + add)) for add in 1:nsamples-1]...,
        SVector{3}(1, nsamples * 2, nsamples + 1),
        [SVector{3}((add + 1, nsamples + add, nsamples + add + 1)) for add in 1:nsamples-1]...,
    ]

    triangles = [t .+ index * nsamples for index in 0:length(dist_eval)-2 for t in single_triangle_set]

    mesh_kwargs = Dict{Symbol, Any}(
        :vertices => vertices,
        :triangles => triangles,
        :components => ones(Int, length(triangles)),
        :number => length(triangles),
    )
    for symbol in Mesh(number=0).unique_keys
        if symbol in keys(mesh_kwargs)
            continue
        end
        mesh_kwargs[symbol] = getproperty(bendy_cylinder, symbol).value
    end

    return Mesh(; mesh_kwargs...)
end

"""
    BendyCylinder(cylinder(s), nsamples=100, height=nothing)

Approximates a cylinder by a [`BendyCylinder`](@ref) with `nsamples` along the circumference.
The `height` of the resulting mesh triangles is by default set to the circumference divided by `nsamples`.
"""
function BendyCylinder(cylinder::Cylinder; nsamples=100, height=nothing)
    if isnothing(height)
        height = 2π * cylinder.radius.value / nsamples
    end
    position = [cylinder.position..., 0.]
    repeats = [cylinder.repeats..., height]

    if isone(cylinder.g_ratio)
        myelin = false
        width = 0.
    else
        myelin = true
        width = cylinder.radius * 2 * (1 - cylinder.g_ratio) / (1 + cylinder.g_ratio)
    end
    iso = cylinder.susceptibility_iso * width
    aniso = cylinder.susceptibility_aniso * width

    last_vec = cross(eachcol(cylinder.rotation)...)
    rotation = [eachcol(cylinder.rotation)..., last_vec]

    bendy_kwargs = Dict{Symbol, Any}(
        :control_point => position,
        :spline_order => 2,
        :repeats => repeats,
        :number => 1,
        :nsamples => nsamples,
        :closed => [0, 0, 1],
        :myelin => myelin,
        :susceptibility_iso => iso,
        :susceptibility_aniso => aniso,
        :rotation => rotation
    )
    for symbol in BendyCylinder(number=0).unique_keys
        if symbol in keys(bendy_kwargs)
            continue
        end
        bendy_kwargs[symbol] = getproperty(cylinder, symbol)
    end
    return BendyCylinder(; bendy_kwargs...)
end

function BendyCylinder(cylinders::Cylinders; nsamples=100, height=nothing)
    if isnothing(height)
        height = 2π * mean(cylinders.radius.value) / nsamples
    end
    BendyCylinder.(cylinders; nsamples=nsamples, height=height)
end


Mesh(c::Cylinder; kwargs...) = Mesh(BendyCylinder(c; kwargs...))
Mesh(c::Cylinders; kwargs...) = Mesh.(BendyCylinder(c; kwargs...))

function Mesh(w::Wall; height=1.)
    normal = w.rotation[:, 1]
    rotation = get_rotation(normal, 3; reference_dimension=:x)
    mesh_kwargs = Dict{Symbol, Any}(
        :vertices => [
            [w.position[1], height/2, height/2],
            [w.position[1], -height/2, height/2],
            [w.position[1], height/2, -height/2],
            [w.position[1], -height/2, -height/2],
        ],
        :triangles => [
            [1, 2, 3],
            [4, 3, 2],
        ],
        :number => 2,
        :rotation => rotation,
        :repeats => [isnothing(w.repeats) ? Inf : w.repeats[1], height, height],
    )
    for symbol in Mesh(number=0).unique_keys
        if symbol in keys(mesh_kwargs)
            continue
        end
        if hasproperty(w, symbol)
            mesh_kwargs[symbol] = getproperty(w, symbol)
        end
    end
    return Mesh(; mesh_kwargs...)
end

Mesh(w::Walls; kwargs...) = Mesh.(w; kwargs...)


"""
    icosahedron(subdivision)

Create an icosahedron.
Based on algorithm from https://danielsieger.com/blog/2021/01/03/generating-platonic-solids.html

Each of the 20 original triangles is sub-divided into subdivision^2 triangles
"""
function icosahedron()
    golden_ratio = (1 + sqrt(5)) / 2
    p = 1 / golden_ratio
    vertices = [v / norm(v) for v in [
        [0, p, -1],
        [p, 1, 0],
        [-p, 1, 0],
        [0, p, 1],
        [0, -p, 1],
        [-1, 0, p],
        [0, -p, -1],
        [1, 0, -p],
        [1, 0, p],
        [-1, 0, -p],
        [p, -1, 0],
        [-p, -1, 0]
    ]]
    triangles = [
        [3, 2, 1],
        [2, 3, 4],
        [6, 5, 4],
        [5, 9, 4],
        [8, 7, 1],
        [7, 10, 1],
        [12, 11, 5],
        [11, 12, 7],
        [10, 6, 3],
        [6, 10, 12],
        [9, 8, 2],
        [8, 9, 11],
        [3, 6, 4],
        [9, 2, 4],
        [10, 3, 1],
        [2, 8, 1],
        [12, 10, 7],
        [8, 11, 7],
        [6, 12, 5],
        [11, 9, 5],
    ]
    return SVector{3, Float64}.(vertices), SVector{3, Int}.(triangles)
end

function icosahedron(subdivisions::Int)
    vertices, base_triangles = icosahedron()
    if isone(subdivisions)
        return vertices, base_triangles
    end

    # place vertices at edges
    edge_map = Dict{Tuple{Int, Int}, Vector{Int}}()
    for triangle in base_triangles
        for (i1, i2) in [
            (triangle[1], triangle[2]),
            (triangle[2], triangle[3]),
            (triangle[3], triangle[1]),
        ]
            edge = (i1, i2)
            if edge in keys(edge_map)
                continue
            end
            edge_map[edge] = length(vertices) .+ (1:(subdivisions - 1))
            edge_map[(i2, i1)] = reverse(edge_map[edge])
            intermediates = [vertices[i1] .* (1 - d) .+ vertices[i2] .* d for d in (0:1/(subdivisions+1):1)[2:end-1]]
            append!(vertices, intermediates)
        end
    end

    triangles = SVector{3, Int}[]
    for triangle in base_triangles
        vertex_indices = Dict{Tuple{Int, Int}, Int}()
        for layer in 1:subdivisions+1
            for index in 1:(subdivisions + 2 - layer)
                if layer == 1
                    if index == 1
                        vertex_indices[(layer, index)] = triangle[1]
                    elseif index == subdivisions + 1
                        vertex_indices[(layer, index)] = triangle[2]
                    else
                        vertex_indices[(layer, index)] = edge_map[(triangle[1], triangle[2])][index - 1]
                    end
                elseif layer == subdivisions + 1
                    vertex_indices[(layer, index)] = triangle[3]
                elseif index == 1
                    vertex_indices[(layer, index)] = edge_map[(triangle[1], triangle[3])][layer - 1]
                elseif layer == subdivisions + 2 - layer
                    vertex_indices[(layer, index)] = edge_map[(triangle[2], triangle[3])][layer - 1]
                else
                    v1 = vertices[edge_map[(triangle[1], triangle[3])][layer - 1]]
                    v2 = vertices[edge_map[(triangle[2], triangle[3])][layer - 1]]
                    dist = (index - 1) / (subdivisions + 1 - layer)
                    push!(vertices, (1 - dist) .* v1 .+ dist .* v2)
                    vertex_indices[(layer, index)] = length(vertices)
                end
            end
        end

        for layer in 1:subdivisions
            for index in 1:(subdivisions + 1 - layer)
                push!(triangles, SVector{3, Int}((
                    vertex_indices[(layer, index)],
                    vertex_indices[(layer, index+1)],
                    vertex_indices[(layer+1, index)],
                )))
            end
        end

        for layer in 1:subdivisions-1
            for index in 2:(subdivisions + 1 - layer)
                push!(triangles, SVector{3, Int}((
                    vertex_indices[(layer, index)],
                    vertex_indices[(layer+1, index)],
                    vertex_indices[(layer+1, index-1)],
                )))
            end
        end
    end

    return [v ./ norm(v) for v in vertices], triangles
end

function Mesh(sphere::Sphere; nsamples=1000)
    subdivisions = Int(ceil(sqrt(nsamples / 20)))
    base_vertices, triangles = icosahedron(subdivisions)

    vertices = [(v .* sphere.radius) .+ sphere.position for v in base_vertices]
    
    mesh_kwargs = Dict{Symbol, Any}(
        :vertices => vertices,
        :triangles => triangles,
        :number => length(triangles),
    )

    for symbol in Mesh(number=0).unique_keys
        if symbol in keys(mesh_kwargs)
            continue
        end
        if hasproperty(sphere, symbol)
            mesh_kwargs[symbol] = getproperty(sphere, symbol)
        end
    end
    return Mesh(; mesh_kwargs...)
end

Mesh(spheres::Spheres; kwargs...) = Mesh.(spheres; kwargs...)

end

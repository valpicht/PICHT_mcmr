module Geometries
using Makie
import StaticArrays: SVector
import LinearAlgebra: cross, ⋅, norm
import Colors
import GeometryBasics
import MCMRSimulator.Plot: PlotPlane, Plot_Geometry, GeometryLike
import MCMRSimulator.Geometries.Internal: FixedGeometry, FixedObstructionGroup, FixedObstruction, Wall, Cylinder, Sphere, obstructions
import MCMRSimulator.Geometries: ObstructionGroup, fix, Mesh, Cylinders, Spheres, Walls


function Makie.plot!(scene::Plot_Geometry{<:Tuple{<:PlotPlane, <:GeometryLike}})
    plot_plane = scene[1]
    base_geometry = scene[2]
    color = scene[:color]
    lc = theme(scene, :linecolor)

    geometry = @lift $base_geometry isa FixedGeometry ? $base_geometry : fix($base_geometry)

    to_plot = @lift project_geometry($plot_plane, $geometry)

    line_color = @lift $color == Makie.automatic ? $lc : $color

    full_kwargs = Dict([key => scene[key] for key in [
        :visible, :overdraw, :fxaa, :transparency, :inspectable, :depth_shift, :model, :space,
        :linestyle, :linewidth
    ]])

    lift(to_plot) do to_iter
        for (func, args, kwargs) in to_iter
            func(scene, args...; color=line_color, kwargs..., full_kwargs...)
        end
    end
    scene
end

Makie.plottype(::PlotPlane, ::GeometryLike) = Plot_Geometry

Makie.convert_arguments(::Plot_Geometry, pp::PlotPlane, geometry::GeometryLike) = (pp, snapshot)


function project_geometry(plot_plane::PlotPlane, geometry::FixedGeometry)
    projections = []
    for t in geometry
        append!(projections, project_geometry(plot_plane, t))
    end
    projections
end

"""
    project_geometry(plot_plane, transform)

Projects the plane on the intrinsic plane of the obstructions deformed by `transform`.
"""
function project_geometry(plot_plane::PlotPlane, group::FixedObstructionGroup{N}) where {N}
    center_obstruction_space = plot_plane.transformation(zero(SVector{3, Float64}))

    obstruction_coordinates_in_plot_plane = SVector{N, SVector{3, Float64}}(map(p->SVector{3, Float64}(plot_plane.transformation(p)) .- center_obstruction_space, eachcol(group.rotation)))

    projections = []
    for obstruction in obstructions(group)
        obstruction_center_in_plot_plane = plot_plane.transformation(group.rotation * obstruction.shift)
        append!(projections, project_obstruction(obstruction.base, obstruction_center_in_plot_plane, obstruction_coordinates_in_plot_plane, group.repeats, (plot_plane.sizex, plot_plane.sizey)))
    end
    projections
end


function project_obstruction(wall::Wall, center::SVector{3, Float64}, obstruction_coordinates_in_plot_plane::SVector{1, SVector{3, Float64}}, repeats::Union{Nothing, SVector{1, Float64}}, sizes::Tuple{<:Real, <:Real})
    constant_center = obstruction_coordinates_in_plot_plane[1] ⋅ center
    normal = obstruction_coordinates_in_plot_plane[1][1:2]
    halfs = (sizes[1]/2, sizes[2]/2)

    function get_line(constant::Float64)
        if abs(normal[1]) > abs(normal[2])
            return [
                SVector{2, Float64}((constant + normal[2] * sizes[2]) / normal[1], -sizes[2]),
                SVector{2, Float64}((constant - normal[2] * sizes[2]) / normal[1], sizes[2]),
                SVector{2, Float64}(NaN, NaN),
            ]
        else
            return [
                SVector{2, Float64}(sizes[1], (constant - normal[1] * sizes[1]) / normal[2]),
                SVector{2, Float64}(-sizes[1], (constant + normal[1] * sizes[1]) / normal[2]),
                SVector{2, Float64}(NaN, NaN),
            ]
        end
    end
    if ~isnothing(repeats)
        repeat_constant_shift = repeats[1] * norm(normal)
        closest_center = mod(constant_center, repeat_constant_shift)
        max_size = abs.(normal) ⋅ halfs
        nshift = div(max_size + abs(closest_center), repeat_constant_shift, RoundUp)
    
        line = SVector{2, Float64}[]
        for shift in -nshift:nshift
            append!(line, get_line(shift * repeat_constant_shift + closest_center))
        end
    else
        line = get_line(constant_center)[1:2]
    end
    [(Makie.lines!, (cut_line(line, halfs), ), Dict())]
end

function project_obstruction(obstruction::Cylinder, center_vec::SVector{3, Float64}, obstruction_coordinates_in_plot_plane::SVector{2, SVector{3, Float64}}, repeats::Union{Nothing, SVector{2, Float64}}, sizes::Tuple{<:Real, <:Real})
    (dirx, diry) = obstruction_coordinates_in_plot_plane
    normal = cross(dirx, diry)
    if normal[3] == 0
        error("Can not plot cylinders or spirals perfectly aligned with plotting plane")
    end
    halfs = (sizes[1]/2, sizes[2]/2)
    center = (center_vec .- normal .* (center_vec[3] / normal[3]))[1:2]
    dirx = (dirx .- normal .* (dirx[3] / normal[3]))[1:2]
    diry = (diry .- normal .* (diry[3] / normal[3]))[1:2]

    projection_size = normal[3] * normal[3]

    theta_normal = atan(normal[2], normal[1])
    function relative_radius(theta)
        ct_rot = cos(theta - theta_normal)
        sqrt(1 / (1 - ct_rot * ct_rot * (1 - projection_size)))
    end

    theta = (0:0.02:2) * π
    radius = obstruction.radius .* relative_radius.(theta)

    points = map((r, t) -> SVector{2, Float64}(r .* cos(t), r * sin(t)), radius, theta)
    push!(points, SVector{2, Float64}(NaN, NaN))

    function get_line(center :: SVector{2, Float64})
        map(p->SVector{2, Float64}(p[1] + center[1], p[2] + center[2]), points)
    end

    if ~isnothing(repeats)
        stepx = dirx * repeats[1]
        stepy = diry * repeats[2]

        toshift = round.([stepx stepy] \ center)
        closest_center = center - (toshift[1] * stepx + toshift[2] * stepy)

        outer_radius = isa(obstruction, Cylinder) ? obstruction.radius : obstruction.outer
        nshiftx = div(abs.(dirx) ⋅ (halfs .+ abs.(closest_center)) + projection_size * outer_radius, repeats[1], RoundUp)
        nshifty = div(abs.(diry) ⋅ (halfs .+ abs.(closest_center)) + projection_size * outer_radius, repeats[2], RoundUp)
        line = SVector{2, Float64}[]
        for i in -nshiftx:nshiftx
            int_pos = closest_center .+ i .* stepx
            for j in -nshifty:nshifty
                pos = SVector{2}(int_pos .+ j .* stepy)
                append!(line, get_line(pos))
            end
        end
    else
        line = get_line(SVector{2}(center))
    end
    [(Makie.lines!, (cut_line(line, halfs), ), Dict())]
end

function cut_line(old_line::AbstractVector{SVector{2, Float64}}, sizes)
    _inside(point::SVector{2, Float64}) = abs(point[1]) <= sizes[1] && abs(point[2]) <= sizes[2]
    new_line = SVector{2, Float64}[]
    prev_point = SVector{2, Float64}(NaN, NaN)

    for new_point in old_line
        if !all(isfinite.(new_point))
            if length(new_line) > 0 && all(isfinite.(new_line[end]))
                push!(new_line, new_point)
            end
        else
            if !all(isfinite.(prev_point))
                if _inside(new_point)
                    push!(new_line, new_point)
                end
            else
                # both prev_point and new_point are finite
                add_overlap!(new_line, prev_point, new_point, sizes)
            end
        end
        prev_point = new_point
    end
    map(p->Makie.Point2(p[1], p[2]), new_line)
end

function add_overlap!(new_line, prev_point, new_point, sizes)
    _inside(point::SVector{2, Float64}) = abs(point[1]) <= sizes[1] && abs(point[2]) <= sizes[2]

    next_inside = _inside(new_point)
    ndiscover = 2 - _inside(prev_point) - next_inside
    if ndiscover == 0
        push!(new_line, new_point)
        return
    end

    # line is at least partially outside
    direction = new_point - prev_point
    dist_new_sq = sum(d->d*d, direction)
    n = norm(direction)
    normal = SVector{2, Float64}(-direction[2]/n, direction[1]/n)
    constant = normal ⋅ prev_point

    discovered = SVector{2, Float64}[]
    function add_valid!(point::SVector{2, Float64})
        distsq = (point - prev_point) ⋅ direction
        if distsq >= 0 && distsq <= dist_new_sq
            push!(discovered, point)
        end
    end
    x1 = (constant - sizes[2] * normal[2]) / normal[1]
    if abs(x1) <= sizes[1]
        add_valid!(SVector{2, Float64}(x1, sizes[2]))
    end
    x2 = (constant + sizes[2] * normal[2]) / normal[1]
    if abs(x2) <= sizes[1]
        add_valid!(SVector{2, Float64}(x2, -sizes[2]))
    end

    y1 = (constant - sizes[1] * normal[1]) / normal[2]
    if abs(y1) < sizes[2]
        add_valid!(SVector{2, Float64}(sizes[1], y1))
    end
    y2 = (constant + sizes[1] * normal[1]) / normal[2]
    if abs(y2) < sizes[2]
        add_valid!(SVector{2, Float64}(-sizes[1], y2))
    end

    if length(discovered) == 0 && ndiscover == 2
        # line is fully outside of the box
        return
    end
    if length(discovered) != ndiscover
        # bad point at edge; let's just ignore it
        return
    end

    if ndiscover == 2
        if length(new_line) > 0 && all(isfinite.(new_line[end]))
            push!(new_line, SVector{2, Float64}(NaN, NaN))
        end
        if (discovered[1] - prev_point) ⋅ direction > (discovered[2] - prev_point) ⋅ direction
            push!(new_line, discovered[2])
            push!(new_line, discovered[1])
        else
            push!(new_line, discovered[1])
            push!(new_line, discovered[2])
        end
    else
        @assert ndiscover == 1
        if next_inside && length(new_line) > 0 && all(isfinite.(new_line[end]))
            push!(new_line, SVector{2, Float64}(NaN, NaN))
        end
        push!(new_line, discovered[1])
        if next_inside
            push!(new_line, new_point)
        end
    end
end

fix_and_mesh(geom::FixedGeometry; height, nsamples) = geom
fix_and_mesh(geom::Mesh; height, nsamples) = fix(geom)
fix_and_mesh(geom::Cylinders; height, nsamples) = fix(Mesh(geom; height=height, nsamples=nsamples == Makie.automatic ? 100 : nsamples))
fix_and_mesh(geom::Spheres; height, nsamples) = fix(Mesh(geom; nsamples=nsamples == Makie.automatic ? 1000 : nsamples))
fix_and_mesh(geom::Walls; height, nsamples) = fix(Mesh(geom; height=height))
fix_and_mesh(geom::ObstructionGroup; height, nsamples) = fix(Mesh(geom))

function Makie.plot!(scene::Plot_Geometry{<:Tuple{<:Nothing, <:GeometryLike}})
    base_geometry = scene[2]
    default_color = scene[:color]
    kwargs = Dict([key => scene[key] for key in [
        :alpha, :visible, :overdraw, :fxaa, :transparency, :inspectable, :depth_shift, :model, :space,
        :shading, :diffuse, :specular, :shininess, :backlight, :ssao
    ]])

    geometry = @lift fix_and_mesh($base_geometry; height=$(scene[:height]), nsamples=$(scene[:nsamples]))

    function plot_group(group, color)
        vert = GeometryBasics.Point{3, Float64}.(group.args.vertices)
        tri = [GeometryBasics.TriangleFace{Int}(o.indices) for o in obstructions(group)]
        patch_color = @lift $default_color == Makie.automatic ? color : $default_color
        geometry_mesh = GeometryBasics.Mesh(vert, tri)
        Makie.mesh!(scene, geometry_mesh; color=patch_color, kwargs...)
    end
    @lift plot_group.($geometry, Colors.distinguishable_colors(length($geometry)))
end

Makie.plottype(::GeometryLike) = Plot_Geometry
end
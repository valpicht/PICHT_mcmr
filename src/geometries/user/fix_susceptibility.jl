"""
Defines [`MCMRSimulator.fix_susceptibility`](@ref MCMRSimulator.Geometries.User.FixSusceptibility.fix_susceptibility).
"""
module FixSusceptibility
import StaticArrays: SVector
import LinearAlgebra: transpose, norm, ⋅, I
import ...Internal.Susceptibility: FixedSusceptibility, SusceptibilityGrid, SusceptibilityGridNoRepeat, SusceptibilityGridRepeat, BaseSusceptibility, CylinderSusceptibility, AnnulusSusceptibility, TriangleSusceptibility, SusceptibilityGridElement, dipole_approximation_repeat, dipole_approximation, IsotropicSusceptibilityGridElement, AnisotropicSusceptibilityGridElement, IsotropicTriangleSusceptibility, AnisotropicTriangleSusceptibility, triangle_magnetisation
import ...Internal.Obstructions.Triangles: FullTriangle, normal, triangle_size
import ...Internal: BoundingBox, radius, lower, upper
import ...Internal.HitGrids: find_hits
import ..Obstructions: ObstructionType, ObstructionGroup, Walls, Cylinders, Spheres, Annuli, Mesh, fields, isglobal, BendyCylinder
import ..SizeScales: grid_resolution

"""
    fix_susceptibility(geometry)

Create a new `FixedSusceptibility` based on the user-defined geometry settings.
"""
function fix_susceptibility(geometry::AbstractVector)
    result = SusceptibilityGrid[]
    for g in geometry
        addition = fix_susceptibility_type(g)
        if ~isnothing(addition)
            push!(result, addition)
        end
    end
    return Tuple(result)
end

fix_susceptibility(group::ObstructionGroup) = fix_susceptibility([group])


fix_susceptibility_type(unknown_group::ObstructionGroup) = nothing

function fix_susceptibility_type(group::Cylinders)
    if all(isone.(group.g_ratio.value))
        return nothing
    end
    b0_field = group.rotation.value[3, :]
    f(a...) = CylinderSusceptibility(a..., b0_field)
    base = f.(group.radius.value, group.g_ratio.value, group.susceptibility_iso.value, group.susceptibility_aniso.value)
    if ~(base isa Vector)
        base = fill(base, length(group))
    end
    add_parent(group, base)
end

function total_susceptibility(group::Cylinders, B0_field::SVector{2, Float64})
    ts = @. group.susceptibility_iso.value + group.susceptibility_aniso.value/4
    gr = group.g_ratio.value
    return @. ts * (1 - gr^2) / (1 + gr)^2 * group.radius.value^2 * 4π
end

function fix_susceptibility_type(group::Annuli)
    if ~any(group.myelin.value)
        return nothing
    end
    b0_field = group.rotation.value[3, :]
    f(a...) = AnnulusSusceptibility(a..., b0_field)
    base = f.(group.inner.value, group.outer.value, group.susceptibility_iso.value, group.susceptibility_aniso.value)
    if ~(base isa Vector)
        base = fill(base, length(group))
    end
    add_parent(group, base; radius_symbol=:outer)
end

function total_susceptibility(group::Annuli, B0_field::SVector{2, Float64})
    ts = @. group.susceptibility_iso.value + group.susceptibility_aniso.value/4
    return @. ts * (group.outer.value^2 - group.inner.value^2) * π
end

for func in (:fix_susceptibility_type, :total_susceptibility)
    @eval $(func)(group::BendyCylinder, args...; kwargs...) = $(func)(Mesh(group), args...; kwargs...)
end

function fix_susceptibility_type(group::Mesh)
    if ~any(group.myelin.value)
        return nothing
    end
    isotropic = all(iszero.(group.susceptibility_aniso.value))
    triangle_type = isotropic ? IsotropicTriangleSusceptibility : AnisotropicTriangleSusceptibility

    b0_field = group.rotation.value[3, :]
    res = triangle_type[]
    radii = Float64[]
    positions = SVector{3, Float64}[]
    for (index, (i1, i2, i3)) in enumerate(group.triangles.value)
        ft = FullTriangle(group.vertices.value[i1], group.vertices.value[i2], group.vertices.value[i3])
        push!(radii, radius(ft) * 4)
        push!(positions, (ft.a + ft.b + ft.c) / 3)
        if isotropic
            push!(res, IsotropicTriangleSusceptibility(ft, group[index].susceptibility_iso))
        else
            push!(res, AnisotropicTriangleSusceptibility(ft, group[index].susceptibility_iso, group[index].susceptibility_aniso, b0_field))
        end
    end
    return add_parent(group, res; positions=positions, radii=radii)
end

function total_susceptibility(mesh::Mesh, B0_field::SVector{3, Float64})
    if all(iszero.(mesh.susceptibility_aniso.value))
        function compute_iso(triangle_index, iso)
            triangle = FullTriangle(mesh.vertices.value[triangle_index]...)
            return iso * triangle_size(triangle)
        end
        return compute_iso.(mesh.triangles.value, mesh.susceptibility_iso.value)
    else
        function compute_aniso(triangle_index, iso, aniso)
            triangle = FullTriangle(mesh.vertices.value[triangle_index]...)
            return triangle_magnetisation(triangle, iso, aniso, B0_field) * triangle_size(triangle)
        end
        return compute_aniso.(mesh.triangles.value, mesh.susceptibility_iso.value, mesh.susceptibility_aniso.value)
    end
end

function add_parent(user::ObstructionGroup, internal::AbstractVector{<:BaseSusceptibility{N}}; positions=nothing, radii=nothing, radius_symbol=:radius) where {N}
    if isnothing(positions)
        positions = isglobal(user.position) ? fill(SVector{N}(user.position.value), length(internal)) : SVector{N}.(user.position.value)
    end

    if isnothing(radii)
        user_radius = getproperty(user, radius_symbol)
        radii = isglobal(user_radius) ? fill(user_radius.value, length(internal)) : user_radius.value
    end

    B0_field = user.rotation.value[3, :]

    susceptibilities = total_susceptibility(user, B0_field)
    if susceptibilities isa Number || susceptibilities isa SVector
        susceptibilities = fill(susceptibilities, length(internal))
    end

    individual_bbs = map(positions, radii) do p, r
        BoundingBox(p .- r, p .+ r)
    end
    if isnothing(user.repeats.value)
        half_repeats = nothing
        bb_indices = BoundingBox(individual_bbs)
    else
        half_repeats = SVector{N}(user.repeats.value) ./ 2
        bb_indices = BoundingBox(-half_repeats, half_repeats)
        bb_off_resonance = bb_indices
    end
    resolution_guess = grid_resolution(user, bb_indices)
    if isinf(resolution_guess)
        orig_size_grid = zero(SVector{N, Int}) .+ 1
    else
        orig_size_grid = Int.(div.(upper(bb_indices) .- lower(bb_indices), resolution_guess, RoundUp))
    end
    resolution = (upper(bb_indices) .- lower(bb_indices)) ./ orig_size_grid
    if isnothing(user.repeats.value)
        size_grid_indices = orig_size_grid .+ 2
        bb_indices = BoundingBox(
            lower(bb_indices) .- resolution,
            upper(bb_indices) .+ resolution,
        )
        size_bb = upper(bb_indices) .- lower(bb_indices)
        size_off_resonance = maximum(size_bb) * 2
        nvoxels_add = @. Int(div((size_off_resonance - size_bb), 2 * resolution, RoundUp))
        bb_off_resonance = BoundingBox(
            lower(bb_indices) .- nvoxels_add .* resolution,
            upper(bb_indices) .+ nvoxels_add .* resolution,
        )
    else
        size_grid_indices = orig_size_grid
    end

    has_hit_bbs = map(positions, radii) do p, r
        half_size = @. max(resolution * 3, r) + 0.5 * resolution
        BoundingBox(p .- half_size, p .+ half_size)
    end

    shifts, grid = find_hits(bb_indices, size_grid_indices, user.repeats.value, has_hit_bbs)

    element_grid = map(grid) do index_arr
        map(index_arr) do (index_obstruction, index_shift)
            (
                SusceptibilityGridElement{N}(
                    positions[index_obstruction],
                    radii[index_obstruction],
                    susceptibilities[index_obstruction]
                ),
                isnothing(user.repeats.value) ? (index=index_obstruction, ) : (
                    index = index_obstruction,
                    shift = index_shift
                ) 
            )
        end
    end

    inv_resolution = 1 ./ resolution

    element_type = eltype(susceptibilities) <: Number ? IsotropicSusceptibilityGridElement{N} : AnisotropicSusceptibilityGridElement{N}

    super_resolution = 3
    if isnothing(user.repeats.value)
        size_grid_off_resonance = @. super_resolution * (size_grid_indices + nvoxels_add * 2)
        grid = zeros(size_grid_off_resonance...)
        @Threads.threads for coordinate in Tuple.(eachindex(IndexCartesian(), grid))
            lower_edge = lower(bb_off_resonance)
            centre = @. ((coordinate - 0.5) * resolution / super_resolution) + lower_edge

            coord_elements = @. div(coordinate - nvoxels_add * super_resolution - 1, super_resolution, RoundDown) + 1
            has_elements = all(coord_elements .>= 1) && all(coord_elements .<= size_grid_indices)
            result = 0.
            for index in 1:length(internal)
                if has_elements && any(e[2].index == index for e in element_grid[coord_elements...])
                    continue
                end
                offset = centre - positions[index]
                result += dipole_approximation(
                    susceptibilities[index],
                    offset,
                    norm(offset),
                    B0_field
                )
            end
            grid[coordinate...] = result
        end
        vec_susceptibilities = eltype(susceptibilities) <: Number ? [s .* B0_field for s in susceptibilities] : susceptibilities
        return SusceptibilityGridNoRepeat{N, eltype(internal), element_type, N*3}(
            inv_resolution,
            transpose(user.rotation.value),
            bb_off_resonance,
            grid,
            super_resolution,
            bb_indices,
            element_grid,
            internal,
            shifts,
            B0_field,
            sum(vec_susceptibilities),
            sum([norm(s) .* p for (s, p) in zip(susceptibilities, positions)]) ./ sum(norm, susceptibilities),
        )
    else
        grid = zeros((super_resolution .* size_grid_indices)...)
        @Threads.threads for coordinate in Tuple.(eachindex(IndexCartesian(), grid))
            centre = @. ((coordinate - 0.5) * resolution / super_resolution) - half_repeats
            result = 0.
            for index in 1:length(internal)
                offset = centre - positions[index]
                result += dipole_approximation_repeat(
                    susceptibilities[index],
                    offset,
                    B0_field,
                    half_repeats .* 2
                )
            end
            coord_elements = @. div(coordinate - 1, super_resolution, RoundDown) + 1
            for (element, source) in element_grid[coord_elements...]
                offset = centre - element.position
                if ~iszero(source.shift)
                    offset = offset .- shifts[source.shift]
                end
                result -= dipole_approximation(
                    element.susceptibility,
                    offset,
                    norm(offset),
                    B0_field,
                )
            end
            grid[coordinate...] = result
        end
        return SusceptibilityGridRepeat{N, eltype(internal), element_type, N*3}(
            inv_resolution,
            transpose(user.rotation.value),
            grid,
            super_resolution,
            half_repeats,
            element_grid,
            internal,
            shifts,
            B0_field,
        )
    end

end

end
"""
Defines [`MCMRSimulator.fix`](@ref MCMRSimulator.Geometries.User.Fix.fix).
"""
module Fix
import LinearAlgebra: I, transpose
import StaticArrays: SVector
import ...Internal: FixedObstructionGroup, FixedGeometry, Internal, get_quadrant, HitGrid, BoundingBox
import ..Obstructions: ObstructionType, ObstructionGroup, Walls, Cylinders, Spheres, Annuli, Mesh, fields, isglobal, BendyCylinder, value_as_vector
import ..SplitMesh: fix_mesh, components
import ..SizeScales: size_scale, grid_resolution

"""
    fix(user_geometry; permeability=0., density=0., dwell_time=0., relaxation=0.)

Creates a fixed version of the user-created geometry that will be used internally by the simulator.
"""
function fix(geometry::AbstractVector; kwargs...)
    result = FixedObstructionGroup[]
    for (original_index, og) in enumerate(geometry)
        append!(result, fix(og, length(result) + 1, original_index; kwargs...))
    end
    return Tuple(result)
end

function fix(geometry::ObstructionGroup, index::Int=1, original_index::Int=1; permeability=0., density=0., dwell_time=0., relaxation=0.)
    for key in propertynames(geometry)
        field_value = getproperty(geometry, key)
        if field_value.field.required && any(isnothing.(field_value.value))
            error("Missing value for required field $field_value")
        end
    end
    result = fix_type(geometry, index, original_index; permeability=permeability, density=density, dwell_time=dwell_time, relaxation=relaxation)
    return result isa FixedObstructionGroup ? (result, ) : Tuple(result)
end


"""
    fix_type(group::ObstructionGroup, index)

Replace user [`ObstructionGroup`](@ref) objects with a [`FixedGeometry`](@ref).

This needs to be defined for each sub-type of [`ObstructionGroup`](@ref),
    e.g. [`Spheres`](@ref), [`Annuli`](@ref), [`Mesh`](@ref).

Shifts, collision properties, and volumetric and surface MRI properties can be added using [`apply_properties`](@ref).
"""
function fix_type end

function fix_type(walls::Walls, index::Int, original_index::Int; kwargs...)
    base_obstructions = fill(Internal.Wall(), length(walls))
    apply_properties(walls, base_obstructions, index, original_index; surface="surface", kwargs...)
end

function fix_type(cylinders::Cylinders, index::Int, original_index::Int; kwargs...)
    if isglobal(cylinders.radius)
        base_obstructions = fill(Internal.Cylinder(cylinders.radius.value), length(cylinders))
    else
        base_obstructions = [Internal.Cylinder(r) for r in cylinders.radius.value]
    end
    apply_properties(cylinders, base_obstructions, index, original_index; surface="surface", volume="inside", kwargs...)
end

function fix_type(spheres::Spheres, index::Int, original_index::Int; kwargs...)
    # Special case for SWC-like substrates represented as overlapping spheres.
    # Instead of converting each sphere independently, we build a single
    # internal object containing all sphere centers and radii.
    if spheres.overlapping.value
        radii = if isglobal(spheres.radius)
            fill(spheres.radius.value, length(spheres))
        else
            spheres.radius.value
        end

        # Overlapping spheres must have explicit positions.
        if !hasproperty(spheres, :position)
            error("Overlapping spheres require positions.")
        end

        positions = if isglobal(spheres.position)
            fill(spheres.position.value, length(spheres))
        else
            spheres.position.value
        end

        # Convert each sphere to (position, radius) format expected internally.
        sphere_data = [
            (SVector{3, Float64}(pos), r)
            for (pos, r) in zip(positions, radii)
        ]

        # Create a single internal geometry object for the whole set
        # of overlapping spheres.
        base_obstructions = [Internal.OverlappingSpheres(sphere_data)]

        # Positions are already encoded in `sphere_data`, so no additional shift
        # should be applied here.
        apply_properties(spheres, base_obstructions, index, original_index;
            surface="surface", volume="inside", apply_shift=false, kwargs...)
    else
        # Standard case: convert each sphere independently.
        if isglobal(spheres.radius)
            base_obstructions = fill(Internal.Sphere(spheres.radius.value), length(spheres))
        else
            base_obstructions = [Internal.Sphere(r) for r in spheres.radius.value]
        end
        apply_properties(spheres, base_obstructions, index, original_index; surface="surface", volume="inside", kwargs...)
    end
end

function fix_type(base_annuli::Annuli, index::Int, original_index::Int; kwargs...)
    annuli = deepcopy(base_annuli)
    annuli.R1_inner_volume = annuli.R1_inner_volume.value .- annuli.R1_outer_volume.value
    annuli.R2_inner_volume = annuli.R2_inner_volume.value .- annuli.R2_outer_volume.value
    annuli.off_resonance_inner_volume = annuli.off_resonance_inner_volume.value .- annuli.off_resonance_outer_volume.value

    if isglobal(annuli.inner)
        base_inner = fill(Internal.Cylinder(annuli.inner.value), length(annuli))
    else
        base_inner = [Internal.Cylinder(r) for r in annuli.inner.value]
    end
    inner_cylinders = apply_properties(annuli, base_inner, index, original_index; surface="inner_surface", volume="inner_volume", kwargs...)

    if isglobal(annuli.outer)
        base_outer = fill(Internal.Cylinder(annuli.outer.value), length(annuli))
    else
        base_outer = [Internal.Cylinder(r) for r in annuli.outer.value]
    end
    outer_cylinders = apply_properties(annuli, base_outer, index + 1, original_index; surface="outer_surface", volume="outer_volume", kwargs...)
    return (inner_cylinders, outer_cylinders)
end


function fix_type(mesh::Mesh, index::Int, original_index::Int; kwargs...)
    return fix_type_single(fix_mesh(mesh), index, original_index; kwargs...)
end

function fix_type_single(mesh::Mesh, index::Int, original_index::Int; kwargs...)
    base_obstructions = [Internal.IndexTriangle(SVector{3, Int32}(index), Int32(c)) for (index, c) in zip(mesh.triangles.value, components(mesh))]
    apply_properties(mesh, base_obstructions, index, original_index; surface="surface", volume="inside", kwargs...)
end

fix_type(bendy_cylinder::BendyCylinder, args...; kwargs...) = fix_type(Mesh(bendy_cylinder), args...; kwargs...)

"""
    apply_properties(user_obstructions, internal_obstructions, index; surface, volume)

Apply the generic obstruction properties defined in `user_obstruction` to `internal_obstructions`.
This function applies (if appropriate):
- positional shifts
- volumetic and surface-bound MRI properties
- collision properties
- rotation
- repeats
- mesh vertices
"""
function apply_properties(user_obstructions::ObstructionGroup, internal_obstructions::Vector{<:Internal.FixedObstruction}, index::Int, original_index::Int; surface=nothing, volume=nothing, apply_shift=true, kwargs...)
    # apply shifts
    if apply_shift && hasproperty(user_obstructions, :position)
        shifts = isglobal(user_obstructions.position) ? fill(user_obstructions.position.value, length(user_obstructions)) : user_obstructions.position.value
        if hasproperty(user_obstructions, :repeats) && ~isnothing(user_obstructions.repeats.value)
            shifts = [mod.(s .+ user_obstructions.repeats.value/2, user_obstructions.repeats.value) .- user_obstructions.repeats.value/2 for s in shifts]
        end
        internal_obstructions = Internal.Shift.(internal_obstructions, shifts)
    end

    fix_array(s::Symbol, other) = other
    fix_array(s::Symbol, ::Nothing) = kwargs[s]
    function fix_array(s::Symbol, arr::Vector{Union{Nothing, T}}) where {T}
        if length(arr) != user_obstructions.n_obstructions
            error("Cannot fix geometry. Parameter $name has $(length(arr)) values, which does not match the number of obstructions ($(user_obstructions.n_obstructions)).")
        end
        [isnothing(v) ? T(kwargs[s]) : v for v in arr]
    end

    # get volumetric MRI properties
    symbols = (:R1, :R2, :off_resonance)
    if ~isnothing(volume)
        get_volume(s) = fix_array(s, getproperty(user_obstructions, Symbol(String(s) * "_" * String(volume))).value)
        values = [get_volume(s) for s in symbols]
        volume = NamedTuple{symbols}(Tuple(values))
    else
        volume = NamedTuple{symbols}(Tuple(zeros(length(symbols))))
    end

    # get surface properties (bound MRI & collision)
    symbols = (:R1, :R2, :off_resonance, :permeability, :density, :dwell_time, :relaxation)
    updated_symbols = replace(symbols, :density=>:surface_density, :relaxation=>:surface_relaxation)
    if ~isnothing(surface)
        get_surface(s) = fix_array(s, getproperty(user_obstructions, Symbol(String(s) * "_" * String(surface))).value)
        values = [get_surface(s) for s in symbols]
        surface = NamedTuple{updated_symbols}(Tuple(values))
    else
        surface = NamedTuple{updated_symbols}(Tuple(zeros(length(symbols))))
    end

    rotation = user_obstructions.rotation.value

    if eltype(internal_obstructions) <: Internal.IndexTriangle
        args = (vertices=SVector{3, Float64}.(user_obstructions.vertices.value), )
    else
        args = NamedTuple()
    end
    bb = BoundingBox(map(o->BoundingBox(o, args), internal_obstructions))
    grid = HitGrid(internal_obstructions, grid_resolution(user_obstructions, bb), user_obstructions.repeats.value, args)

    if eltype(internal_obstructions) <: Internal.IndexTriangle
        args = (
            inside_mask=BitArray(undef, (maximum(components(user_obstructions)), size(grid.indices)...)), 
            inside_prepared=Ref(false),
            args...
        )
    end
    
    result = FixedObstructionGroup(
        internal_obstructions,
        user_obstructions.repeats.value,
        index,
        original_index,
        rotation,
        grid,
        volume,
        surface,
        size_scale(user_obstructions),
        args
    )
    return result
end


end
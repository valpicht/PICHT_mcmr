"""
Defines the properties of a certain geometry type.

Types:
- `ObstructionType`

Functions:
- `fields`
- `key_value_pairs`
"""
module ObstructionTypes

import StaticArrays: MVector, SMatrix
import LinearAlgebra: I
import ..Fields: Field, property_fields, FieldValue

"""
Defines individual obstruction types
"""
struct ObstructionType{N}
    singular :: Symbol
    plural :: Symbol
    ndim :: Int
    surfaces :: Vector{Symbol}
    volumes :: Vector{Symbol}
    group_volumes :: Bool
    extra_fields :: Vector{Field}
    include_shift :: Bool
end

function ObstructionType(singular::Symbol; plural=nothing, ndim=3, surfaces=[:surface], volumes=[:inside], group_volumes=false, fields=Field[], include_shift=true)
    if isnothing(plural)
        plural = Symbol(String(singular) * "s")
    end
    ObstructionType{plural}(singular, plural, ndim, surfaces, volumes, group_volumes, fields, include_shift)
end

function Base.show(io::IO, field::ObstructionType)
    print(io, String(field.singular))
end

function fields(ot::ObstructionType)
    position_type = ot.ndim == 1 ? Float64 : MVector{ot.ndim, Float64}
    default_rotation = SMatrix{3, ot.ndim}(I(3)[:, 1:ot.ndim])
    fields = [
        ot.extra_fields...,
        Field{SMatrix{3, ot.ndim, Float64, 3 * ot.ndim}}(:rotation, "Rotation applied to all obstructions in group. Can be set to a matrix or one of :x, :y, or, :z (see [`MCMRSimulator.get_rotation`](@ref MCMRSimulator.Methods.get_rotation)).", default_rotation, only_group=true, required=true),
        Field{MVector{ot.ndim, Float64}}(:repeats, "Length scale on which the obstructions are repeated (um).", only_group=true),
        Field{Bool}(:use_boundingbox, "Use bounding boxes for an initial filtering of possible intersections.", true, only_group=true),
        Field{Float64}(:grid_resolution, "Resolution of the grid that the volume is split up into (um). Defaults to roughly one grid element per obstruction.", only_group=true),
        property_fields...,
    ]
    if ot.include_shift
        push!(fields, Field{position_type}(:position, "Spatial offset of obstruction from origin.", zero(position_type), required=true))
    end

    fields = filter(f -> (
        ~(f.per_surface || f.per_volume) ||
        ~iszero(length(categories(f, ot)))
        ), fields)
    sort!(fields, by=f -> ~f.required)
    return fields
end

function categories(field::Field, ot::ObstructionType)
    categories = []
    if field.per_surface
        append!(categories, ot.surfaces)
    end
    if field.per_volume
        append!(categories, ot.volumes)
    end
    return categories
end


"""
    key_value_pairs(obstruction_type, n_obstructions)

Defines a mapping from symbols to [`FieldValue`](@ref) objects.

Returns a type with:
- a sub-set of symbols that are enough to set all unique field values.
- a mapping from symbols to field values.
"""
function key_value_pairs(field::Field, ot::ObstructionType, n_obstructions=0)
    D = Dict{Symbol, FieldValue}
    if ~field.per_surface && ~field.per_volume
        return [field.name], D(field.name => FieldValue(field, n_obstructions))
    end
    cs = categories(field, ot)
    values = [FieldValue(field, n_obstructions, c) for c in cs]
    long_name(ct) = Symbol(String(field.name) * "_" * String(ct))
    long_name2(ct) = Symbol(String(ct) * "_" * String(field.name))
    if length(cs) == 0
        return Symbol[], D()
    elseif length(cs) == 1
        return [field.name], D(field.name => values[1], long_name(cs[1]) => values[1], long_name2(cs[1]) => values[1])
    else
        return long_name.(cs), D(
            [long_name(c) => v for (c, v) in zip(cs, values)]...,
            [long_name2(c) => v for (c, v) in zip(cs, values)]...
        )
    end
end

function key_value_pairs(fields::Vector{<:Field}, ot::ObstructionType, n_obstructions=0)
    (unique_keys, mapping) = zip([
        key_value_pairs(f, ot, n_obstructions) for f in fields
    ]...)

    return vcat(unique_keys...), mergewith((a, b) -> error("Duplicate field name found"), mapping...)
end

function key_value_pairs(ot::ObstructionType, n_obstructions=0)
    return key_value_pairs(fields(ot), ot, n_obstructions)
end

end
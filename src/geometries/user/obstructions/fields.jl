"""
Defines a flexible `Field` interface to set the properties of user-defined geometries.
"""
module Fields

import StaticArrays: SMatrix
import .....Methods: get_rotation

"""
    Field{T}(name::Symbol[, description::String[, default_value::T]]; per_surface=false, per_volume=false, only_group=false, required=false)

Represents a property of an obstruction with:
- `Type`: type of the property value (can also be Nothing unless `required` is true).
- `name`: name of the property.
- `description`: describes the property to the user.
- `default_value`: what the value is set to if not provided by the user.
- `per_surface`: can have a different value for each surface in the obstruction.
- `per_volume`: can have a different value within each enclosed volume of the obstruction.
- `required`: if true, the user has to provide a value, which cannot be nothing.
"""
struct Field{T}
    name :: Symbol
    description :: String
    default_value :: Union{Nothing, T}
    per_surface :: Bool
    per_volume :: Bool
    only_group :: Bool
    required :: Bool
    function Field{T}(name, description="", default_value=nothing; per_surface=false, per_volume=false, only_group=false, required=false) where {T}
        new{T}(name, description, default_value, per_surface, per_volume, only_group, required)
    end
end

description(f::Field) = f.description

"""
    convert_value(T, value, required::Bool)

Convert `value` into type `T` or a vector of `T`.
If required is false `value` can also be nothing
or a vector containing `T` and nothing.
"""
function convert_value(::Type{T}, value) where {T}
    if isnothing(value) || value isa T
        return value
    end
    try
        return T(value)
    catch
        if ~(value isa AbstractVector)
            error("New values can not be converted into $T or vector of $T's: $value")
        end
        if value isa Vector{Union{Nothing, T}}
            return value
        elseif value isa AbstractVector{<:Union{Nothing, T}}
            return Vector{Union{Nothing, T}}(value)
        else
            return Union{Nothing, T}[isnothing(v) ? v : T(v) for v in value]
        end
    end
end

convert_value(T::Type{SMatrix{3, N, Float64, M}}, value) where {N, M} = T(get_rotation(value, N))

"""
    FieldValue(field::Field{T}, n_obstructions, [, category::Symbol][, value])

Represents the value for a specific field across a vector of obstructions.
For fields with multiple values (one per volume or surface),
`category` is required to set this.
The current value is given by `value`.
This has to be type `T` or a vector of type `T`.
For a non-required field, it can also be nothing (or a vector containing nothings).
"""
mutable struct FieldValue{T}
    field :: Field
    n_obstructions :: Int
    category :: Union{Nothing, Symbol}
    value :: Union{Nothing, T, Vector{Union{Nothing, T}}}
    function FieldValue(field::Field{T}, n_obstructions, varargs...) where {T}
        if field.per_surface || field.per_volume
            (category, varargs...) = varargs
        else
            category = nothing
        end
        if iszero(length(varargs))
            value = nothing
        elseif isone(length(varargs))
            value = varargs[1]
        else
            error("Too many arguments to FieldValue")
        end
        if field.per_surface || field.per_volume
            @assert category isa Symbol
        else
            @assert isnothing(category)
        end
        if isnothing(value)
            value = field.default_value
        end
        res = new{T}(field, n_obstructions, category, value)
        return res
    end
end

value_as_vector(fv::FieldValue{T}) where {T} = fv.value isa Union{Nothing, T} ? fill(fv.value, fv.n_obstructions) : fv.value

description(fv::FieldValue{T}) where {T} = (
    description(fv.field) * 
    (isnothing(fv.category) ? "" : " $(uppercasefirst(String(fv.category))) property.") *
    (fv.field.required ? " Field is required." : " Field can be null.") *
    (" Expected type: $T.")
)

function Base.setproperty!(v::FieldValue{T}, s::Symbol, value) where {T}
    if s == :value
        value = convert_value(T, value)
        if v.field.required && any(isnothing.(value))
            error("Cannot set nothing(s) to required $(v.field)")
        end
        if v.field.only_group && ~isglobal(v, value)
            error("Cannot set an array of values to $(v.field), which expects to be the same across obstructions.")
        end
        if v.n_obstructions > 0 && ~isglobal(v, value) && length(value) != v.n_obstructions
            error("Expected one value for each of the $(v.n_obstructions) obstructions to $(v.field), not $(length(value)) values.")
        end
    end
    setfield!(v, s, value)
end

isglobal(fv::FieldValue{T}) where {T} = fv.value isa Union{Nothing, T}
isglobal(::FieldValue{T}, value) where {T} = value isa Union{Nothing, T}

Base.length(fv::FieldValue) = fv.n_obstructions
Base.getindex(fv::FieldValue{T}, index::Integer) where {T} = isglobal(fv) ? fv.value : fv.value[index]
function Base.setindex!(fv::FieldValue{T}, value, index::Integer) where {T}
    if fv.field.only_group
        error("$(fv.field) only accepts a single value across a whole obstruction group.")
    end
    if isglobal(fv)
        fv.value = fill(fv.value, fv.n_obstructions)
    end
    if isnothing(value)
        if fv.field.required
            error("Cannot set nothing to obstruction in required $(fv.field)")
        end
    elseif ~(value isa T)
        value = T(value)
    end
    fv.value[index] = value
end


function Base.show(io::IO, field::Field)
    print(io, "Field($(field.name))")
end

function Base.show(io::IO, fv::FieldValue)
    text = isnothing(fv.category) ? "" : " in $(fv.category)"
    print(io, "$(fv.field)$(text) = $(fv.value)")
end

property_fields = (
    Field{Float64}(:R1, "Additional longitudinal relaxation rate (kHz).", 0., per_volume=true, per_surface=true, required=true),
    Field{Float64}(:R2, "Additional transverse relaxation rate (kHz).", 0., per_volume=true, per_surface=true, required=true),
    Field{Float64}(:off_resonance, "Additional off-resonance field offset (kHz).", 0., per_volume=true, per_surface=true, required=true),
    Field{Float64}(:dwell_time, "Average time a particle stays stuck to the surface (ms).", per_surface=true),
    Field{Float64}(:density, "Surface density of stuck particles relative to the volume density (um).", per_surface=true),
    Field{Float64}(:permeability, "Rate of particle passing through the obstruction in arbitrary units.", per_surface=true),
    Field{Float64}(:relaxation, "Rate of signal loss at each collision. The actual signal loss at each collision is e^(-x * sqrt(t)), where x is this rate and t is the timestep.", per_surface=true),
    Field{Float64}(:size_scale, "Size of the smallest obstructions. If not set explicitly, this will be determined by the minimum radius or distance between objects (see `size_scale`).", only_group=true),
)

end
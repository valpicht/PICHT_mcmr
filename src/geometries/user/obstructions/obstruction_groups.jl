"""
Defines `ObstructionGroup`, which is the parent type of all user-defined geometries.
"""
module ObstructionGroups

import ..ObstructionTypes: ObstructionType, fields, key_value_pairs
import ..Fields: Field, FieldValue, isglobal


"""
Parent type of all geometry objects like [`Cylinders`](@ref MCMRSimulator.Geometries.User.Obstructions.Cylinders) or [`Mesh`](@ref MCMRSimulator.Geometries.User.Obstructions.Mesh).
"""
struct ObstructionGroup{N}
    type :: ObstructionType{N}
    field_values :: Dict
    unique_keys :: Vector{Symbol}
    n_obstructions :: Int
end


Base.length(g::ObstructionGroup) = g.n_obstructions
function Base.getindex(g::ObstructionGroup, index::Int)
    return IndexedObstruction(g, index)
end
Base.iterate(g::ObstructionGroup) = iterate(g, 1)
Base.iterate(g::ObstructionGroup, state::Int) = state > g.n_obstructions ? nothing : (g[state], state + 1)

"""
A single instance of an [`ObstructionGroup`](@ref) obtained by indexing.
"""
struct IndexedObstruction{N}
    group :: ObstructionGroup{N}
    index :: Int
end


Base.propertynames(og::ObstructionGroup) = og.unique_keys
Base.propertynames(io::IndexedObstruction) = propertynames(io.group)
function Base.getproperty(og::ObstructionGroup, s::Symbol) 
    d = getfield(og, :field_values)
    if s in keys(d)
        return d[s]
    else
        return getfield(og, s)
    end
end

function Base.setproperty!(og::ObstructionGroup, s::Symbol, value) 
    og.field_values[s].value = value
end

function Base.getproperty(io::IndexedObstruction, s::Symbol) 
    if s in (:group, :index)
        return getfield(io, s)
    end
    return io.group.field_values[s][io.index]
end

function Base.setproperty!(io::IndexedObstruction, s::Symbol, value) 
    io.group.field_values[s][io.index] = value
end

"""
    nvolumes(group)

For any `ObstructionType` with group_volumes, this computes the number of distinct volumes.

This is mainly used to compute the number of distinct components in a mesh.
"""
function nvolumes(group::ObstructionGroup)
    if length(group.type.volumes) == 0
        error("Getting the number of volumes is not defiend for $(group.type)")
    end
    if group.type.group_volumes
        return 1
    else
        return group.n_obstructions
    end
end

function ObstructionGroup(type::ObstructionType; number=nothing, kwargs...)

    (unique_keys, field_values) = key_value_pairs(type, isnothing(number) ? 0 : number)
    for (key, value) in kwargs
        field_values[key].value = value
    end

    if isnothing(number)
        # try to guess the number of obstructions
        for fv in values(field_values)
            if fv.field.per_volume && type.group_volumes
                continue
            end
            if ~isglobal(fv)
                if !isnothing(number) && length(fv.value) != number
                    error("Inconsistent number of obstructions inferred when constructing $(type). Are there $number or $(length(fv.value)) elements?")
                end
                number = length(fv.value)
            end
        end
        if isnothing(number)
            number = 1
        end
    end

    group = ObstructionGroup(type, field_values, unique_keys, number)
    for fv in values(field_values)
        if fv.field.per_volume && type.group_volumes
            fv.n_obstructions = nvolumes(group)
        else
            fv.n_obstructions = number
        end
        if ~isglobal(fv) && (length(fv.value) != fv.n_obstructions)
            error("Incorrect number of values for $fv. Expected $(fv.n_obstructions); got $(length(fv.value)).")
        end
    end
    return group
end

function Base.show(io::IO, og::ObstructionGroup)
    print(io, "$(og.n_obstructions) $(lowercase(String(og.type.plural)))(")
    not_set = Symbol[]
    required_missing = Symbol[]
    for key in propertynames(og)
        fv = og.field_values[key]
        if isnothing(fv.value)
            if fv.field.required
                push!(required_missing, key)
            else
                push!(not_set, key)
            end
        elseif (isglobal(fv) && key != :vertices) || og.n_obstructions < 5
            print(io, "$key = $(fv.value), ")
        else
            print(io, "$key = <unique values>, ")
        end
    end
    print(io, ")")
    if length(required_missing) > 0
        text = join(required_missing, ", ", " and ")
        print(io, "\nMissing required fields: $text")
    end
    if length(not_set) > 0
        text = join(not_set, ", ", " and ")
        print(io, "\nNot set: $text")
    end
end

function Base.show(io::IO, obstruction::IndexedObstruction)
    og = obstruction.group
    print(io, "$(og.type.singular) at index $(obstruction.index) (")
    already_done = Set()
    for (key, fv) in og.field_values
        if fv in already_done
            continue
        end
        push!(already_done, fv)
        my_value = fv[obstruction.index]
        if ~isnothing(my_value)
            print(io, "$key = $(my_value), ")
        end
    end
    print(io, ")")
end

end
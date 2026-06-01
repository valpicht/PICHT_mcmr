"""
Supports I/O to and from json

Functions:
- `write_geometry`
- `read_geometry`
"""
module JSON
import JSON as JSON_pkg
import ..Obstructions: Obstructions, ObstructionGroup, description

"""
    write_geometry([io, ]geometry)

Writes the user-defined geometry as a JSON file.
If no IO is provided, the geometry is written to stdout.
"""
write_geometry(object) = write_geometry(stdout, object)

function write_geometry(io::IO, geometry::AbstractVector{<:ObstructionGroup}) 
    print(io, "[\n")
    for group in geometry
        write_geometry(io, group)
        if group != last(geometry)
            print(io, ",")
        end
        print(io, "\n")
    end
    print(io, "]\n")
end

function write_geometry(io::IO, group::ObstructionGroup)
    print(io, "  {\n")
    print(io, "    \"type\": \"$(group.type.plural)\",\n")
    print(io, "    \"number\": $(group.n_obstructions),\n")
    for unique_key in group.unique_keys
        fv = group.field_values[unique_key]
        print(io, "    \"#$(unique_key)_description\": \"$(description(fv))\",\n")
        print(io, "    \"$unique_key\": ")
        print(io, JSON_pkg.json(fv.value))
        if unique_key != last(group.unique_keys)
            print(io, ",")
        end
        print(io, "\n")
    end
    print(io, "  }\n")
end

function write_geometry(filename::AbstractString, geometry)
    open(io -> write_geometry(io, geometry), filename; write=true)
end

"""
    read_geometry(filename/io)

Read geometry from a JSON file.
"""
function read_geometry(io::IO)
    obj = JSON_pkg.parse(io)
    parse_geometry(obj)
end

function read_geometry(filename::AbstractString)
    if startswith(filename, "  {\n") || startswith(filename, "[\n")
        obj = JSON_pkg.parse(filename)
        return parse_geometry(obj)
    end
    open(read_geometry, filename; read=true)
end

function parse_geometry(v::Vector)
    return parse_geometry.(v)
end

function parse_geometry(d::Dict)
    type_name = Symbol(d["type"])
    constructor = getproperty(Obstructions, type_name)
    kwargs = Dict{Symbol, Any}()
    for (key, value) in d
        if key == "type"
            continue
        end
        if startswith(key, "#") && endswith(key, "_description")
            continue
        end
        kwargs[Symbol(key)] = value
    end
    return constructor(; kwargs...)
end


end
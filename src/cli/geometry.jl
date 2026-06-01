"""
Defines command line interface for `mcmr geometry`
"""
module Geometry

import ArgParse: ArgParse, ArgParseSettings, @add_arg_table!, add_arg_table!, parse_args, ArgParseError, usage_string
import StaticArrays: MMatrix, MVector
import Random
import CSV
import Tables
import ...Geometries.User.Obstructions: Walls, Spheres, Cylinders, Annuli, BendyCylinder, fields, field_to_docs, Field, ObstructionGroup, FieldValue
import ...Geometries.User.JSON: write_geometry, read_geometry
import ...Geometries.User.RandomDistribution: random_positions_radii
import ...Methods: get_rotation


field_type(::Field{T}) where {T} = T

struct FieldParser{T}
    value :: Union{T, Vector{T}}
end

value_as_vector(fp::FieldParser{T}) where {T} = fp.value isa T ? [fp.value] : fp.value

function ArgParse.parse_item(::Type{FieldParser{T}}, text::AbstractString) where {T<:Number}
    try
        return FieldParser{T}(parse(T, text))
    catch
    end
    for dlm in (',', ':')
        if dlm in text
            try
                return FieldParser{T}([parse(T, x) for x in split(text, dlm)])
            catch
            end
        end
    end
    as_mat = CSV.read(text, Tables.matrix, delim=' ', ignorerepeated=true, header=false)
    if size(as_mat) == (1, 1)
        return FieldParser{T}(T(as_mat[1, 1]))
    elseif size(as_mat, 1) == 1
        as_vec = as_mat[1, :]
    elseif size(as_mat, 2) == 1
        as_vec = as_mat[:, 1]
    else
        error("Expected just a row or column of numbers in $(text), got a whole table instead.")
    end
    return FieldParser{T}(T.(as_vec))
end

function ArgParse.parse_item(::Type{FieldParser{MVector{1, T}}}, text::AbstractString) where {T}
    sub_parse = ArgParse.parse_item(FieldParser{T}, text).value
    if sub_parse isa T
        FieldParser{MVector{1, T}}(MVector{1, T}([sub_parse]))
    else
        FieldParser{MVector{1, T}}([MVector{1, T}([s]) for s in sub_parse])
    end
end

function ArgParse.parse_item(::Type{FieldParser{T}}, text::AbstractString) where {T<:AbstractVector}
    if ':' in text && ',' in text
        return FieldParser{T}([T([parse(eltype(T), x) for x in split(sub, ',')]) for sub in split(text, ':')])
    end
    for dlm in (',', ':')
        if dlm in text
            return FieldParser{T}(T([parse(eltype(T), x) for x in split(text, dlm)]))
        end
    end
    as_mat = CSV.read(text, Tables.matrix, delim=' ', ignorerepeated=true, header=false)
    return FieldParser{T}(T.(eachrow(as_mat)))
end

struct RotationParser{N}
    value :: MMatrix{3, N, Float64}
end

function ArgParse.parse_item(::Type{RotationParser{N}}, text::AbstractString) where {N}
    if length(text) == 1
        base = Symbol(text)
    else
        try
            if ',' in text
                base = parse.(Float64, split(text, ','))
            end
        catch
            as_mat = CSV.read(text, Tables.matrix, delim=' ', ignorerepeated=true, header=false)
            if size(as_mat, 1) == 1
                base = as_mat[1, :]
            elseif size(as_mat, 2) == 1
                base = as_mat[:, 1]
            else
                base = as_mat
            end
        end
    end
    return RotationParser(get_rotation(base, N))
end

"""
    get_parser()

Returns the parser of arguments for `mcmr geometry`
"""
function get_parser(; kwargs...)
    parser = ArgParseSettings(
        prog="mcmr geometry", 
        description="Various commands to manipulate the geometry used in MCMR simulations";
        kwargs...
    )
    @add_arg_table! parser begin
        "create"
            help="Creates a geometry JSON file containing a single set of obstructions (e.g., cylinders, spheres, walls)."
            action=:command
        "create-random"
            help="Creates a geometry JSON file containing randomly distributed cylinders or spheres."
            action=:command
        "merge"
            help="Merge the geometries in multiple JSON files into one."
            action=:command
    end

    parser["create"].description = "Create a geometry JSON file with obstructions of a specific type. Select a type for more help."

    for (as_string, constructor) in [
        ("walls", Walls),
        ("cylinders", Cylinders),
        ("spheres", Spheres),
        ("annuli", Annuli),
        ("bendy-cylinder", BendyCylinder),
    ]
        for sub_command in ("create", "create-random")
            if sub_command == "create-random" && as_string == "walls"
                continue
            end
            add_arg_table!(parser[sub_command], as_string, Dict(
                :help => "Fill the geometry with $as_string.",
                :action => :command
            ))
            parser[sub_command][as_string].description = "Create a geometry JSON file filled with only $as_string with any properties defined by the flags."


            if sub_command == "create"
                add_arg_table!(parser[sub_command][as_string],
                    "number", Dict(
                        :arg_type => Int,
                        :help => "Number of obstructions to create.",
                        :required => true,
                    )
                )
            elseif sub_command == "create-random"
                add_arg_table!(parser[sub_command][as_string],
                    "target-density", Dict(
                        :arg_type => Float64,
                        :help => "Total fraction of the total space that should be filled with these obstructions.",
                        :required => true,
                    ),
                    "--mean-radius", Dict(
                        :arg_type => Float64,
                        :help => "Mean radius of the $as_string (um).",
                        :default => 1.,
                    ),
                    "--var-radius", Dict(
                        :arg_type => Float64,
                        :help => "Variance of the $as_string radius distribution (um^2).",
                        :default => 0.,
                    ),
                    "--seed", Dict(
                        :arg_type => Int64,
                        :help => "Initialisation for random number seed. Supply this to get reproducible results.",
                    ),
                )
                if as_string == "annuli"
                    add_arg_table!(parser[sub_command][as_string],
                        "--g-ratio", Dict(
                            :arg_type => Float64,
                            :help => "Ratio of the inner over outer radius for the annuli.",
                            :default => 0.7,
                        )
                    )
                end
            end
            add_arg_table!(parser[sub_command][as_string],
                "output_file", Dict(
                    :help => "Geometry JSON output filename.",
                    :required => true,
                )
            )
            group = constructor(number=0)
            for unique_key in group.unique_keys
                field_value = group.field_values[unique_key]
                as_dict = Dict(
                    :dest_name => String(unique_key),
                    :help => field_to_docs(field_value),
                    :required => field_value.field.required && isnothing(field_value.field.default_value),
                    :arg_type => FieldParser{field_type(field_value.field)},
                )
                if ~isnothing(field_value.field.default_value)
                    as_dict[:default] = FieldParser{field_type(field_value.field)}(field_value.field.default_value)
                end
                flag = "--" * String(unique_key)
                if field_type(field_value.field) == Bool
                    for s in (:required, :default, :arg_type)
                        pop!(as_dict, s)
                    end
                    if field_value.field.default_value
                        flag = "--no-" * String(unique_key)
                        as_dict[:action] = :store_false
                    else
                        as_dict[:action] = :store_true
                    end
                elseif field_value.field.name == :rotation
                    as_dict[:arg_type] = RotationParser{group.type.ndim}
                    as_dict[:default] = RotationParser{group.type.ndim}(get_rotation(:I, group.type.ndim))
                end

                if sub_command == "create-random"
                    if unique_key in [:radius, :position, :inner, :outer]
                        continue
                    end
                    if unique_key == :repeats
                        as_dict[:required] = true
                        as_dict[:help] = "Length scale on which the random pattern of $as_string repeats itself (um)."
                    end
                end
                add_arg_table!(parser[sub_command][as_string], flag, as_dict)
            end
        end
    end

    add_arg_table!(parser["merge"],
        "output_file", Dict(
            :help => "A new geometry JSON file containing all the obstructions from the input files.",
            :required => true,
        ),
        "input_file", Dict(
            :help => "List of input geometry JSON files that should be merged. Each input file can contain obstructions of different types.",
            :required => true,
            :nargs => '+',
        ),
    )
    return parser
end

"""
    parse_user_argument(value, n_objects)

Parse the user argument to something that can actually be used internally by Julia.
"""
parse_user_argument(value, n_objects) = value
function parse_user_argument(value::FieldParser{T}, n_objects) where {T}
    if value.value isa T
        return value.value
    else
        @assert length(value.value) == n_objects
        return value.value
    end
end
parse_user_argument(value::RotationParser, n_objects) = value.value


function parse_user_rotation(value::Vector)
    if length(value) == 1
        return Symbol(value[1])
    end
    parsed = parse.(Float64, value)
    if length(parsed) == 3
        return parsed
    else
        return reshape(parsed, (length(parsed) // 3, 3))
    end
end

"""
    run_main([arguments])

Runs the `mcmr geometry` command line interface.
The supplied arguments can be provided as a sequence of strings (as provided by the terminal)
or as a dictionary (as provided by `ArgParse` after parsing).
By default it is set to `ARGS`.
"""
function run_main(args=ARGS::AbstractVector[<:AbstractString]; kwargs...)
    parser = get_parser(; kwargs...)
    return run_main(parse_args(args, parser))
end

function run_main(::Nothing)
    return Cint(1)
end

function run_main(args::Dict{<:AbstractString, <:Any})
    cmd = args["%COMMAND%"]
    if cmd == "create"
        run_create(args[cmd])
    elseif cmd == "create-random"
        run_create_random(args[cmd])
    elseif cmd == "merge"
        run_merge(args[cmd])
    end
    return Cint(0)
end

function run_create(args::Dict{<:AbstractString, <:Any})
    obstruction_type = args["%COMMAND%"]
    flags = args[obstruction_type]
    output_file = pop!(flags, "output_file")
    constructor = Dict(
        "walls" => Walls,
        "cylinders" => Cylinders,
        "spheres" => Spheres,
        "annuli" => Annuli,
        "bendy-cylinder" => BendyCylinder,
    )[obstruction_type]
    number = pop!(flags, "number")
    symbol_flags = Dict(Symbol(k) => parse_user_argument(v, number) for (k, v) in flags)
    filtered = Dict(k=>v for (k, v) in symbol_flags if ~isnothing(v))
    result = constructor(;number=number, filtered...)
    write_geometry(output_file, result)
end

function run_create_random(args::Dict{<:AbstractString, <:Any})
    obstruction_type = args["%COMMAND%"]
    flags = args[obstruction_type]
    if "seed" in keys(flags)
        Random.seed!(pop!(flags, "seed"))
    end
    output_file = pop!(flags, "output_file")
    (constructor, ndim) = Dict(
        "cylinders" => (Cylinders, 2),
        "spheres" => (Spheres, 3),
        "annuli" => (Annuli, 2),
    )[obstruction_type]
    (positions, radius) = random_positions_radii(
        flags["repeats"].value, pop!(flags, "target-density"), ndim; 
        mean=pop!(flags, "mean-radius"), variance=pop!(flags, "var-radius"), min_radius=0.
    )
    flags["position"] = positions
    if obstruction_type == "annuli"
        flags["outer"] = radius
        flags["inner"] = radius * pop!(flags, "g-ratio")
    else
        flags["radius"] = radius
    end
    number = length(radius)
    symbol_flags = Dict(Symbol(k) => parse_user_argument(v, number) for (k, v) in flags)
    filtered = Dict(k=>v for (k, v) in symbol_flags if ~isnothing(v))
    result = constructor(;number=number, filtered...)
    write_geometry(output_file, result)
end

function run_merge(args)
    res = ObstructionGroup[]
    for input_file in args["input_file"]
        input_geom = read_geometry(input_file)
        if input_geom isa ObstructionGroup
            push!(res, input_geom)
        else
            append!(res, input_geom)
        end
    end
    write_geometry(args["output_file"], res)
end

end
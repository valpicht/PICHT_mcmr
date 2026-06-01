"""
Defines the command line interface to `mcmr run`.
"""
module Run

import ArgParse: ArgParseSettings, @add_arg_table!, add_arg_group!, parse_args
import DataFrames: DataFrame
import CSV
import Tables
import Random
import MRIBuilder: read_sequence, Sequence
import ...Geometries.User.JSON: read_geometry
import ...Simulations: Simulation
import ...Spins: Snapshot, BoundingBox, longitudinal, transverse, phase, orientation, position
import ...Evolve: readout
import ...Subsets: Subset, get_subset


"""Add simulation definition parameters to an argument parser."""
function add_simulation_definition!(parser)
    add_arg_group!(parser, "Define the simulation parameters", :simulation)
    @add_arg_table! parser begin
        "geometry"
            required = true
            help = "JSON file describing the spatial configuration of any obstructions as well as biophysical properties associated with those obstructions. Can be generated using `mcmr geometry`. Alternatively, a mesh file can be provided."
        "sequence"
            nargs = '*'
            help = "One of more pulseq .seq files describing the sequences to be run."
        "-D", "--diffusivity"
            help = "Diffusivity of free water (um^2/ms)."
            arg_type = Float64
            default = 3.
        "--R1"
            help = "Longitudinal relaxation in 1/ms. This relaxation rate will at the very least be applied to free, extra-cellular spins. It might be overriden in the 'geometry' for bound spins or spins inside any obstructions."
            arg_type = Float64
            default = 0.
        "--R2"
            help = "Transverse relaxation in 1/ms. This relaxation rate will at the very least be applied to free, extra-cellular spins. It might be overriden in the 'geometry' for bound spins or spins inside any obstructions."
            arg_type = Float64
            default = 0.
        "--bvecs"
            help = "ASCII text file with gradient orientations in FSL format (https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FDT/UserGuide#Processing_pipeline)."
    end
end

"""Add output flags to an argument parser."""
function add_output_flags!(parser)
    add_arg_group!(parser, "Output flags. At least one is required", :output, required=true)
    @add_arg_table! parser begin
        "--output-signal", "-o"
            help = "Writes the total signal at the readouts to this file as a comma-separated value (CSV) table."
        "--output-snapshot"
            help = "Writes the state of all the spins at the readouts to this file as a comma-separated value (CSV) table."
    end
end

"""Add readout flags to an argument parser."""
function add_readout_flags!(parser)
    add_arg_group!(parser, "Readout flags. These control when the signal/spin states will be read out", :readout)
    @add_arg_table! parser begin
        "--nTR"
            help = "Acquire the signal provided at the sequence readouts for this many repetition times (TRs). Output will be stored as a CSV file."
            arg_type = Int
            default = 1
        "--times", "-T"
            help = "Acquire the signal at the given times within each TR (in ms). Multiple values can be provided (e.g., '-T 0 10 15.3'). By default, the Readout markers in the sequence will be used instead."
            arg_type = Float64
            nargs = '+'
        "--skip-TR"
            help = "The number of repetition times the simulation will run before starting to acquire data."
            arg_type = Int
            default = 0
        "--subset"
            help = """Can be provided multiple times. For each time it is provided, the signal will be computed at each readout for a specific subset of spins.  This subset is defined by one or two values from bound/free/inside/outside. Afterwards they can include an integer value to select a specific geometry to consider the bound/inside state of. An additional integer value could be given to select a specific obstruciton within that geometry.
            For example:
            - `--subset free`: include any free spins
            - `--subset inside`: include any spins inside any geometry
            - `--subset outside 2`: include any spins outside of the second obstruction group in the geometry
            - `--subset inside bound 2 3`: include any spins stuck to the inside surface of the 3rd obstruction in the second obstruction group of the geometry.
            """
            nargs = '+'
            action = :append_arg
    end
end

"""Add initialisation flags to an argument parser."""
function add_init_flags!(parser)
    add_arg_group!(parser, "Initialisation flags. These control the spins initial state", :init)
    @add_arg_table! parser begin
        "--Nspins", "-N"
            help = "Number of spins to simulate. Ignored if --init is set."
            arg_type = Int
            default = 10000
        "--voxel-size"
            help = "Size of the voxel (in mm) over which the initial spins are spread."
            arg_type = Float64
            default = 1.
        "--longitudinal"
            help = "Initial value of the longitudinal magnetisation for each spin. Note the the equilibrium longitudinal magnetisation for each spin is 1."
            arg_type = Float64
            default = 1.
        "--transverse"
            help = "Initial value of the magnitude of the transverse magnetisation for each spin."
            arg_type = Float64
            default = 0.
        "--phase"
            help = "Initial value of the phase of the transverse magnetisation for each spin in degrees."
            arg_type = Float64
            default = 0.
        #"--init"
        #    help = "Continues the simulation with the spin states stored in a JSON file (produced using the 'mcmr pre-run' command). If used all other initialisation flags (except for --seed) are ignored."
        "--seed"
            help = "Initialisation for random number seed. Supply this to get reproducible results. If --init is also set, this flag will override the seed stored in this initialisation file."
            arg_type = Int
    end
end

"""
    get_parser()

Returns the parser of arguments for `mcmr run`
"""
function get_parser(; kwargs...)
    parser = ArgParseSettings(prog="mcmr run", description="Runs a Monte Carlo simulation of the MRI signal evolution for spins interacting with the geometry."; kwargs...)
    add_simulation_definition!(parser)
    add_output_flags!(parser)
    add_readout_flags!(parser)
    add_init_flags!(parser)
    return parser
end


"""
    read_bvecs(filename)

Reads the bvecs into a Nx3 matrix.
"""
function read_bvecs(filename)
    table = CSV.read(filename, Tables.matrix, delim=' ', ignorerepeated=true, header=false)
    (r, c) = size(table)
    if c != 3
        if r != 3
            error("bvec shape is $(size(table)) rather than (N, 3) or (3, N)")
        end
        table = transpose(table)
    end
    return table
end

"""
    run_main([arguments])

Runs the `mcmr run` command line interface.
The supplied arguments can be provided as a sequence of strings (as provided by the terminal)
or as a dictionary (as provided by `ArgParse` after parsing).
By default it is set to `ARGS`.
"""
function run_main(args=ARGS::AbstractVector[<:AbstractString]; kwargs...)
    parser = get_parser(; kwargs...)
    run_main(parse_args(args, parser))
end

run_main(::Nothing) = Cint(1)

function run_main(args::Dict{<:AbstractString, <:Any})
    if "seed" in keys(args)
        Random.seed!(args["seed"])
    end
    geometry = read_geometry(args["geometry"])
    sequences = read_sequence.(args["sequence"])

    if isnothing(args["bvecs"])
        sequence_indices = 1:length(sequences)
        bvec_indices = zeros(Int, length(sequences))
        all_sequences = sequences
    else
        sequence_indices = Int[]
        bvec_indices = Int[]
        all_sequences = Sequence[]
        bvecs = read_bvecs(args["bvecs"])
        for (index_seq, sequence) in enumerate(sequences)
            if can_rotate_bvec(sequence)
                append!(sequence_indices, fill(index_seq, size(bvecs, 1)))
                append!(bvec_indices, 1:size(bvecs, 1))
                append!(all_sequences, [rotate_bvec(sequence, bv) for bv in eachrow(bvecs)])
            else
                push!(sequence_indices, index_seq)
                push!(bvec_indices, 0)
                push!(all_sequences, sequence)
            end
        end
        if all(iszero.(bvec_indices))
            @warn "None of the input sequences include bvec-dependent gradients, so the `--bvec` flag will have no effect."
        end
    end

    sequence_names = [args["sequence"][i] for i in sequence_indices]

    simulation = Simulation(all_sequences; geometry=geometry, R1=args["R1"], R2=args["R2"], diffusivity=args["diffusivity"])

    # initial state
    bb = BoundingBox(args["voxel-size"] * 1000/2)
    init_snapshot = Snapshot(args["Nspins"], simulation, bb; longitudinal=args["longitudinal"], transverse=args["transverse"])
    as_snapshot = !isnothing(args["output-snapshot"])
    readout_times = iszero(length(args["times"])) ? nothing : args["times"]
    subsets = [Subset(), parse_subset.(args["subset"])...]
    result = readout(init_snapshot, simulation, readout_times; skip_TR=args["skip-TR"], nTR=args["nTR"], noflatten=true, return_snapshot=as_snapshot, subset=subsets)

    # convert to tabular format
    if !isnothing(args["output-signal"])
        df_list = []
        for index in eachindex(IndexCartesian(), result)
            if as_snapshot
                value = SpinOrientationSum(result[index])
            else
                value = result[index]
            end
            orient_as_vec = orientation(value)
            push!(df_list, (
                sequence=sequence_names[index[1]],
                sequence_index=sequence_indices[index[1]],
                bvec=bvec_indices[index[1]],
                TR=index[3],
                readout=index[2],
                subset=index[4] - 1,
                nspins=length(value),
                longitudinal=longitudinal(value),
                transverse=transverse(value),
                phase=phase(value),
                Sx=orient_as_vec[1],
                Sy=orient_as_vec[2],
            ))
        end
        df = DataFrame(df_list)
        @show df
        CSV.write(args["output-signal"], df)
    end

    if !isnothing(args["output-snapshot"])
        df_list = []
        for index in eachindex(IndexCartesian(), result)
            snapshot = result[index]
            for (ispin, spin) in enumerate(snapshot)
                orient_as_vec = orientation(spin)
                pos =  position(spin)
                push!(df_list, (
                    sequence=sequence_names[index[1]],
                    sequence_index=sequence_indices[index[1]],
                    bvec=bvec_indices[index[1]],
                    TR=index[3],
                    readout=index[2],
                    spin=ispin,
                    x=pos[1],
                    y=pos[2],
                    z=pos[3],
                    longitudinal=longitudinal(spin),
                    transverse=transverse(spin),
                    phase=phase(spin),
                    Sx=orient_as_vec[1],
                    Sy=orient_as_vec[2],
                ))
            end
        end
        df = DataFrame(df_list)
        CSV.write(args["output-snapshot"], df)
    end

    return Cint(0)
end


function parse_subset(arguments)
    arg, others... = arguments
    kwargs = Dict{Symbol, Any}()
    if arg == "bound"
        kwargs[:bound] = true
    elseif arg == "free"
        kwargs[:bound] = false
    elseif arg == "inside"
        kwargs[:inside] = true
    elseif arg == "outside"
        kwargs[:inside] = false
    else
        error("First argument of --subset should be bound, free, inside, or outside; not $inside_str")
    end
    if length(others) == 0
        return Subset(; kwargs...)
    end

    arg, others... = others
    if arg in ["bound", "free"]
        if :bound in keys(kwargs)
            error("Second argument of --subset cannot be bound or free, if the first argument is already bound or free.")
        end
        kwargs[:bound] = arg == "bound"
    elseif arg in ["inside", "outside"]
        if :inside in keys(kwargs)
            error("Second argument of --subset cannot be inside or outside, if the first argument is already inside or outside.")
        end
        kwargs[:inside] = arg == "inside"
    else
        try
            kwargs[:geometry_index] = parse(Int, arg)
        catch
            error("Second argument should be one of bound/free/inside/outside or parse to an integer. Instead, '$index' was received.")
        end
        if length(others) > 1
            error("Too many arguments provided to --subset.")
        elseif length(others) == 1
            kwargs[:obstruction_index] = parse(Int, others[1])
        end
        return Subset(kwargs...)
    end

    if length(others) == 0
        return Subset(kwargs...)
    end
    arg, others... = others
    kwargs[:geometry_index] = parse(arg, Int)

    if length(others) == 0
        return Subset(kwargs...)
    end
    arg, others... = others
    kwargs[:obstruction_index] = parse(arg, Int)
    if length(others) > 0
        error("Too many arguments provided to --subset.")
    end
    @assert length(others) == 0
    return Subset(kwargs...)
end

end
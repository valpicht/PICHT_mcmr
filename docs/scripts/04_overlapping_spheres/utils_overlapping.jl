# ============================================================
# utils_overlapping.jl
#
# Utility functions for overlapping-sphere substrates.
#
# This file is included by the scripts in:
#
#     docs/scripts/04_overlapping_spheres/
#
# The goal is to keep the main scripts short and readable by
# collecting repeated helper functions here.
# ============================================================

using MCMRSimulator
using MRIBuilder
using DelimitedFiles

# ============================================================
# File paths
# ============================================================

"""
    default_substrate_dir()

Return the directory where example substrate files are expected to be stored.

By default, SWC files and example meshes are stored in:

    substrate_generation/example_outputs/

This avoids using absolute paths such as:

    /home/valentine/project/Grow_neurons/results/...
"""
function default_substrate_dir()
    return joinpath(
        "substrate_generation",
        "example_outputs",
    )
end

"""
    add_swc_extension_if_missing(filename::String)

Return a valid `.swc` filename.

If the input string already ends with `.swc`, it is returned unchanged.
Otherwise, the `.swc` extension is appended automatically.
"""
function add_swc_extension_if_missing(filename::String)
    if endswith(filename, ".swc")
        return filename
    else
        return filename * ".swc"
    end
end

"""
    swc_path_from_name(swc_name::String)

Construct the full relative path to an SWC file stored in the default
substrate directory.

The input can be given either with or without the `.swc` extension.

Example:

    swc_path_from_name("one_soma_list_new")

returns:

    "substrate_generation/example_outputs/one_soma_list_new.swc"
"""
function swc_path_from_name(swc_name::String)
    filename = add_swc_extension_if_missing(swc_name)

    return joinpath(
        default_substrate_dir(),
        filename,
    )
end

"""
    check_file_exists(path::String; description = "file")

Check that a file exists and throw a clear error message if it does not.
"""
function check_file_exists(path::String; description = "file")
    if !isfile(path)
        error("""
        Missing $description.

        Expected path:
            $path

        Please copy the required file into:
            $(default_substrate_dir())
        """)
    end

    return true
end

# ============================================================
# SWC loading
# ============================================================

"""
    read_swc_spheres(swc_path::String)

Read an SWC file and extract sphere positions and radii.

Expected SWC format:

    id type x y z radius parent

Lines starting with `#` are ignored.

Returns:
- `positions`: vector of 3-element position vectors `[x, y, z]`;
- `radii`: vector of sphere radii.

This is the representation needed to create an overlapping-sphere
substrate in MCMR.
"""
function read_swc_spheres(swc_path::String)
    check_file_exists(swc_path; description = "SWC file")

    positions = Vector{Vector{Float64}}()
    radii = Float64[]

    open(swc_path, "r") do io
        for line in eachline(io)
            stripped = strip(line)

            if isempty(stripped) || startswith(stripped, "#")
                continue
            end

            fields = split(stripped)

            if length(fields) < 7
                @warn "Skipping malformed SWC line" line
                continue
            end

            x = parse(Float64, fields[3])
            y = parse(Float64, fields[4])
            z = parse(Float64, fields[5])
            radius = parse(Float64, fields[6])

            push!(positions, [x, y, z])
            push!(radii, radius)
        end
    end

    if isempty(positions)
        error("""
        No spheres were read from the SWC file.

        File:
            $swc_path

        Please check that the file uses the expected SWC format:
            id type x y z radius parent
        """)
    end

    println("SWC file loaded successfully.")
    println("Number of spheres: ", length(radii))
    println("Minimum radius: ", minimum(radii))
    println("Maximum radius: ", maximum(radii))

    return positions, radii
end

"""
    load_overlapping_spheres(swc_name::String)

Load an SWC file from the example substrate directory and convert it into
an MCMR `Spheres` geometry with `overlapping = true`.

Example:

    geometry, positions, radii = load_overlapping_spheres("one_soma_list_new")

The `.swc` extension is optional.
"""
function load_overlapping_spheres(swc_name::String)
    swc_path = swc_path_from_name(swc_name)

    println("SWC name: ", swc_name)
    println("SWC path: ", swc_path)

    positions, radii = read_swc_spheres(swc_path)

    geometry = Spheres(
        radius = radii,
        position = positions,
        overlapping = true,
    )

    println("Overlapping-sphere geometry created.")
    println("Geometry type: ", typeof(geometry))

    return geometry, positions, radii
end

# ============================================================
# Geometry inspection and bounding box
# ============================================================

"""
    print_sphere_ranges(positions, radii)

Print the x, y, and z ranges of a set of spheres.

The ranges include the sphere radii, so the printed extent corresponds to
the full substrate extent, not only the sphere centers.
"""
function print_sphere_ranges(positions, radii)
    xs = [p[1] for p in positions]
    ys = [p[2] for p in positions]
    zs = [p[3] for p in positions]

    xmin = minimum(xs .- radii)
    xmax = maximum(xs .+ radii)

    ymin = minimum(ys .- radii)
    ymax = maximum(ys .+ radii)

    zmin = minimum(zs .- radii)
    zmax = maximum(zs .+ radii)

    println("Substrate spatial extent:")
    println("x range: ", xmin, " to ", xmax)
    println("y range: ", ymin, " to ", ymax)
    println("z range: ", zmin, " to ", zmax)
end

"""
    make_bbox_from_spheres(positions, radii; margin = 2.0)

Create a bounding box around a set of spheres.

The bounding box includes the radius of each sphere and adds an optional
margin around the complete substrate.
"""
function make_bbox_from_spheres(positions, radii; margin = 2.0)
    xs = [p[1] for p in positions]
    ys = [p[2] for p in positions]
    zs = [p[3] for p in positions]

    xmin = minimum(xs .- radii)
    xmax = maximum(xs .+ radii)

    ymin = minimum(ys .- radii)
    ymax = maximum(ys .+ radii)

    zmin = minimum(zs .- radii)
    zmax = maximum(zs .+ radii)

    return BoundingBox(
        [xmin - margin, ymin - margin, zmin - margin],
        [xmax + margin, ymax + margin, zmax + margin],
    )
end

"""
    create_spins_for_spheres(positions, radii; nspins = 10_000, margin = 2.0)

Create a random spin snapshot in a bounding box around the sphere substrate.
"""
function create_spins_for_spheres(positions, radii; nspins = 10_000, margin = 2.0)
    bbox = make_bbox_from_spheres(
        positions,
        radii;
        margin = margin,
    )

    println("Bounding box:")
    println(bbox)

    spins = Snapshot(nspins, bbox)

    println("Number of generated spins: ", length(spins))

    return spins
end

"""
    print_inside_outside_counts(geometry, spins)

Print how many spins are inside and outside the geometry.

This is useful as a sanity check before running a simulation.
"""
function print_inside_outside_counts(geometry, spins)
    inside_mask = isinside(geometry, spins) .> 0

    n_inside = count(inside_mask)
    n_outside = length(spins) - n_inside

    println("Inside spins: ", n_inside)
    println("Outside spins: ", n_outside)

    return n_inside, n_outside
end

# ============================================================
# MRI sequence and simulation
# ============================================================

"""
    create_dwi_sequence(; bval = 2.0, TE = 80, TR = 300)

Create a diffusion-weighted imaging sequence using `MRIBuilder`.
"""
function create_dwi_sequence(; bval = 2.0, TE = 80, TR = 300)
    return DWI(
        bval = bval,
        TE = TE,
        TR = TR,
        scanner = MRIBuilder.Siemens_Prisma,
    )
end

"""
    create_dwi_simulation(geometry; bval = 2.0, TE = 80, TR = 300, diffusivity = 1.0)

Create a diffusion-weighted simulation on an overlapping-sphere geometry.
"""
function create_dwi_simulation(
    geometry;
    bval = 2.0,
    TE = 80,
    TR = 300,
    diffusivity = 1.0,
)
    sequence = create_dwi_sequence(
        bval = bval,
        TE = TE,
        TR = TR,
    )

    return Simulation(
        sequence,
        diffusivity = diffusivity,
        geometry = geometry,
    )
end

"""
    readout_inside_outside(spins, simulation; skip_TR = 2)

Run a readout and split the result into inside and outside compartments.

Returns:
- `inside_signal`
- `outside_signal`
"""
function readout_inside_outside(spins, simulation; skip_TR = 2)
    signals = readout(
        spins,
        simulation,
        skip_TR = skip_TR,
        subset = [
            Subset(inside = true),
            Subset(inside = false),
        ],
    )

    inside_signal = signals[1]
    outside_signal = signals[2]

    return inside_signal, outside_signal
end

"""
    signal_scalar(signal)

Extract a scalar signal value from an MCMR signal object.

The function returns the magnitude of the transverse signal.
"""
function signal_scalar(signal)
    return abs(signal.orient.transverse)
end

"""
    simulate_signal_vs_b(spins, geometry, bvals; TE = 80, TR = 300, diffusivity = 1.0, skip_TR = 2)

Run one DWI simulation per b-value and return the inside and outside
scalar transverse signals.

Returns:
- `inside_signals`
- `outside_signals`
"""
function simulate_signal_vs_b(
    spins,
    geometry,
    bvals;
    TE = 80,
    TR = 300,
    diffusivity = 1.0,
    skip_TR = 2,
)
    inside_signals = Float64[]
    outside_signals = Float64[]

    for b in bvals
        println("Running simulation for TE = $TE ms, b = $b")

        simulation = create_dwi_simulation(
            geometry;
            bval = b,
            TE = TE,
            TR = TR,
            diffusivity = diffusivity,
        )

        inside_signal, outside_signal = readout_inside_outside(
            spins,
            simulation;
            skip_TR = skip_TR,
        )

        push!(inside_signals, signal_scalar(inside_signal))
        push!(outside_signals, signal_scalar(outside_signal))
    end

    return inside_signals, outside_signals
end

# ============================================================
# Signal processing
# ============================================================

"""
    normalize_by_first_value(signals)

Normalize a signal vector by its first value, usually S(0).
"""
function normalize_by_first_value(signals)
    s0 = signals[1]

    if isnan(s0) || s0 == 0
        return fill(NaN, length(signals))
    end

    return signals ./ s0
end

"""
    log_normalized_signal(signals)

Compute ln(S/S0), safely handling NaN or non-positive values.
"""
function log_normalized_signal(signals)
    s0 = signals[1]
    ln_norm = Float64[]

    for s in signals
        if !isnan(s) && !isnan(s0) && s > 0 && s0 > 0
            push!(ln_norm, log(s / s0))
        else
            push!(ln_norm, NaN)
        end
    end

    return ln_norm
end

# ============================================================
# CSV output
# ============================================================

"""
    ensure_output_dirs()

Create the standard output directories used by the scripts.
"""
function ensure_output_dirs()
    output_dir = joinpath("docs", "data", "processed")
    figure_dir = joinpath("docs", "figures")

    mkpath(output_dir)
    mkpath(figure_dir)

    return output_dir, figure_dir
end

"""
    write_signal_vs_b_csv(output_file, bvals, inside_signals, outside_signals)

Save raw and normalized inside/outside signals for a b-value sweep.
"""
function write_signal_vs_b_csv(
    output_file,
    bvals,
    inside_signals,
    outside_signals,
)
    inside_norm = normalize_by_first_value(inside_signals)
    outside_norm = normalize_by_first_value(outside_signals)

    data = hcat(
        bvals,
        inside_signals,
        inside_norm,
        outside_signals,
        outside_norm,
    )

    open(output_file, "w") do io
        writedlm(
            io,
            ["b" "inside_signal" "inside_norm" "outside_signal" "outside_norm"],
            ',',
        )
        writedlm(io, data, ',')
    end

    println("Signal data saved:")
    println(output_file)
end
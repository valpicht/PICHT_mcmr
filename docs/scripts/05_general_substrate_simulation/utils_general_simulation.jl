# ============================================================
# utils_general_simulation.jl
#
# Utility functions for general single-substrate DWI simulations.
#
# This file is included by scripts in:
#
#     docs/scripts/05_general_substrate_simulation/
#
# The goal is to provide reusable functions for loading substrates
# from SWC or CSV files, converting them into overlapping-sphere
# geometries, running DWI simulations, computing ln(S/S0), estimating
# ADC, and saving results.
# ============================================================

using MCMRSimulator
using MRIBuilder
using DelimitedFiles
using Statistics
using Random

# ============================================================
# Paths and file handling
# ============================================================

"""
    default_substrate_dir()

Return the folder where example substrate files are stored.

Expected location:

    substrate_generation/example_outputs/

This avoids hard-coded absolute paths such as:

    /home/valentine/project/Grow_neurons/results/...
"""
function default_substrate_dir()
    return joinpath(
        "substrate_generation",
        "example_outputs",
    )
end

"""
    ensure_output_dirs()

Create and return the standard output directories used by the scripts.

Returns:
- `output_dir`: `docs/data/processed`
- `figure_dir`: `docs/figures`
"""
function ensure_output_dirs()
    output_dir = joinpath("docs", "data", "processed")
    figure_dir = joinpath("docs", "figures")

    mkpath(output_dir)
    mkpath(figure_dir)

    return output_dir, figure_dir
end

"""
    add_extension_if_missing(filename::String, extension::String)

Return `filename` with `extension` appended if missing.

Examples:

    add_extension_if_missing("substrate", ".swc")
    add_extension_if_missing("substrate.swc", ".swc")
"""
function add_extension_if_missing(filename::String, extension::String)
    if endswith(filename, extension)
        return filename
    else
        return filename * extension
    end
end

"""
    substrate_path_from_name(substrate_name::String, substrate_format::String)

Return the path to a substrate file in the default substrate directory.

Supported formats:
- `"swc"`
- `"csv"`

Examples:

    substrate_path_from_name("one_soma_list_new", "swc")
    substrate_path_from_name("my_substrate", "csv")
"""
function substrate_path_from_name(substrate_name::String, substrate_format::String)
    format = lowercase(substrate_format)

    if format == "swc"
        filename = add_extension_if_missing(substrate_name, ".swc")
    elseif format == "csv"
        filename = add_extension_if_missing(substrate_name, ".csv")
    else
        error("""
        Unsupported substrate format: $substrate_format

        Supported formats are:
            "swc"
            "csv"
        """)
    end

    return joinpath(
        default_substrate_dir(),
        filename,
    )
end

"""
    check_file_exists(path::String; description = "file")

Check that a file exists and throw a clear error if not.
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

"""
    print_selected_substrate(substrate_name, substrate_format)

Print information about the selected substrate.
"""
function print_selected_substrate(substrate_name::String, substrate_format::String)
    path = substrate_path_from_name(substrate_name, substrate_format)

    println("Selected substrate")
    println("------------------")
    println("name:   ", substrate_name)
    println("format: ", substrate_format)
    println("path:   ", path)
end

# ============================================================
# Substrate loading: SWC and CSV
# ============================================================

"""
    read_swc_spheres(swc_path::String)

Read an SWC file and extract sphere positions and radii.

Expected SWC format:

    id type x y z radius parent

Lines starting with `#` are ignored.

Returns:
- `positions`: vector of 3-element vectors `[x, y, z]`
- `radii`: vector of sphere radii
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

        Expected SWC format:
            id type x y z radius parent
        """)
    end

    return positions, radii
end

"""
    read_csv_spheres(csv_path::String)

Read a CSV file and extract sphere positions and radii.

Expected columns:

    x, y, z, radius

This function is intentionally permissive:
- it ignores an optional header if the first row contains text;
- it accepts comma-separated numeric rows.

Returns:
- `positions`: vector of 3-element vectors `[x, y, z]`
- `radii`: vector of sphere radii
"""
function read_csv_spheres(csv_path::String)
    check_file_exists(csv_path; description = "CSV file")

    positions = Vector{Vector{Float64}}()
    radii = Float64[]

    open(csv_path, "r") do io
        for (line_number, line) in enumerate(eachline(io))
            stripped = strip(line)

            if isempty(stripped) || startswith(stripped, "#")
                continue
            end

            fields = split(stripped, ",")

            if length(fields) < 4
                @warn "Skipping malformed CSV line" line_number line
                continue
            end

            try
                x = parse(Float64, strip(fields[1]))
                y = parse(Float64, strip(fields[2]))
                z = parse(Float64, strip(fields[3]))
                radius = parse(Float64, strip(fields[4]))

                push!(positions, [x, y, z])
                push!(radii, radius)

            catch
                if line_number == 1
                    # Most likely a header row.
                    continue
                else
                    @warn "Skipping non-numeric CSV line" line_number line
                end
            end
        end
    end

    if isempty(positions)
        error("""
        No spheres were read from the CSV file.

        File:
            $csv_path

        Expected CSV columns:
            x,y,z,radius
        """)
    end

    return positions, radii
end

"""
    load_substrate_spheres(substrate_name::String, substrate_format::String)

Load sphere positions and radii from a substrate file.

Supported formats:
- `"swc"`
- `"csv"`

Returns:
- `positions`
- `radii`
- `substrate_path`
"""
function load_substrate_spheres(substrate_name::String, substrate_format::String)
    path = substrate_path_from_name(substrate_name, substrate_format)
    format = lowercase(substrate_format)

    println("Loading substrate from:")
    println("  ", path)

    if format == "swc"
        positions, radii = read_swc_spheres(path)
    elseif format == "csv"
        positions, radii = read_csv_spheres(path)
    else
        error("Unsupported substrate format: $substrate_format")
    end

    println("Substrate loaded successfully.")
    println("Number of spheres: ", length(radii))
    println("Minimum radius: ", minimum(radii))
    println("Maximum radius: ", maximum(radii))

    return positions, radii, path
end

"""
    create_overlapping_spheres_geometry(positions, radii)

Create an MCMR `Spheres` geometry with `overlapping = true`.
"""
function create_overlapping_spheres_geometry(positions, radii)
    geometry = Spheres(
        radius = radii,
        position = positions,
        overlapping = true,
    )

    println("Overlapping-sphere geometry created.")
    println("Geometry type: ", typeof(geometry))

    return geometry
end

"""
    load_overlapping_substrate(substrate_name::String, substrate_format::String)

Load a substrate from SWC or CSV and convert it into an overlapping-sphere
MCMR geometry.

Returns:
- `geometry`
- `positions`
- `radii`
- `substrate_path`
"""
function load_overlapping_substrate(substrate_name::String, substrate_format::String)
    positions, radii, substrate_path = load_substrate_spheres(
        substrate_name,
        substrate_format,
    )

    geometry = create_overlapping_spheres_geometry(
        positions,
        radii,
    )

    return geometry, positions, radii, substrate_path
end

# ============================================================
# Geometry inspection and spin generation
# ============================================================

"""
    print_sphere_ranges(positions, radii)

Print the x, y, and z ranges of a set of spheres, including their radii.
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
    make_bbox_from_positions(positions, radii; margin = 2.0)

Create a bounding box around a set of spheres.

The bounding box includes the radius of each sphere and adds an optional
margin around the full substrate.
"""
function make_bbox_from_positions(positions, radii; margin = 2.0)
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
    create_spins(positions, radii; nspins = 10_000, margin = 2.0, seed = 1234)

Create random spins in a bounding box around the substrate.

Returns:
- `spins`
- `bbox`
"""
function create_spins(
    positions,
    radii;
    nspins = 10_000,
    margin = 2.0,
    seed = 1234,
)
    Random.seed!(seed)

    bbox = make_bbox_from_positions(
        positions,
        radii;
        margin = margin,
    )

    println("Bounding box:")
    println(bbox)

    spins = Snapshot(nspins, bbox)

    println("Number of generated spins: ", length(spins))

    return spins, bbox
end

"""
    inside_outside_counts(geometry, spins)

Count how many spins are inside and outside a geometry.

Returns:
- `n_inside`
- `n_outside`
"""
function inside_outside_counts(geometry, spins)
    inside_mask = isinside(geometry, spins) .> 0

    n_inside = count(inside_mask)
    n_outside = length(spins) - n_inside

    return n_inside, n_outside
end

"""
    print_inside_outside_counts(geometry, spins)

Print and return the number of inside and outside spins.
"""
function print_inside_outside_counts(geometry, spins)
    n_inside, n_outside = inside_outside_counts(geometry, spins)

    println("Inside spins: ", n_inside)
    println("Outside spins: ", n_outside)

    return n_inside, n_outside
end

"""
    split_inside_outside_spins(geometry, spins)

Split a spin snapshot into inside and outside spins.

Returns:
- `inside_spins`
- `outside_spins`
"""
function split_inside_outside_spins(geometry, spins)
    inside_mask = isinside(geometry, spins) .> 0

    inside_spins = spins[inside_mask]
    outside_spins = spins[.!inside_mask]

    return inside_spins, outside_spins
end

# ============================================================
# DWI sequences and simulations
# ============================================================

"""
    get_directional_sequences(bvals, TE, direction; TR = 300)

Create a set of DWI sequences for one diffusion direction.

A reference sequence is created at the maximum b-value. Lower b-values
are obtained by rescaling the diffusion gradients using `MRIBuilder.adjust`.

This follows the workflow used in the original simulation scripts.
"""
function get_directional_sequences(bvals, TE, direction; TR = 300)
    if maximum(bvals) == 0
        error("maximum(bvals) must be larger than 0.")
    end

    ref_sequence = DWI(
        bval = maximum(bvals),
        TE = TE,
        TR = TR,
        gradient = (
            orientation = direction,
        ),
        scanner = MRIBuilder.Siemens_Prisma,
    )

    return MRIBuilder.adjust(
        ref_sequence,
        diffusion = (
            scale = sqrt.(bvals ./ maximum(bvals)),
        ),
        merge = false,
    )
end

"""
    run_directional_readout(spins, geometry, bvals, TE, direction; TR, diffusivity, timestep, skip_TR)

Run a combined DWI simulation for one TE and one diffusion direction.

Returns:
- vector of scalar transverse signal values, one per b-value.
"""
function run_directional_readout(
    spins,
    geometry,
    bvals,
    TE,
    direction;
    TR = 300,
    diffusivity = 1.0,
    timestep = nothing,
    skip_TR = 0,
)
    sequences = get_directional_sequences(
        bvals,
        TE,
        direction;
        TR = TR,
    )

    if isnothing(timestep)
        simulation = Simulation(
            sequences,
            diffusivity = diffusivity,
            geometry = geometry,
        )
    else
        simulation = Simulation(
            sequences,
            diffusivity = diffusivity,
            geometry = geometry,
            timestep = timestep,
        )
    end

    sigs = readout(
        spins,
        simulation,
        skip_TR = skip_TR,
    )

    signals = [abs(sig.orient.transverse) for sig in sigs]

    return signals
end

"""
    run_single_TE_direction(geometry, spins, bvals, TE, direction_name, direction; kwargs...)

Run one simulation block for one TE and one diffusion direction.

Returns:
- `signals`
- `ln_signal`
- `adc`
- `intercept`
- `fit_ln`
"""
function run_single_TE_direction(
    geometry,
    spins,
    bvals,
    TE,
    direction_name,
    direction;
    TR = 300,
    diffusivity = 1.0,
    timestep = nothing,
    skip_TR = 0,
)
    println("\nRunning TE = $TE ms, direction = $direction_name")
    println("Gradient orientation = ", direction)

    signals = @time run_directional_readout(
        spins,
        geometry,
        bvals,
        TE,
        direction;
        TR = TR,
        diffusivity = diffusivity,
        timestep = timestep,
        skip_TR = skip_TR,
    )

    ln_signal = log_normalized_signal(signals)
    adc, intercept = fit_adc_from_ln_signal(bvals, ln_signal)
    fit_ln = fitted_ln_signal(bvals, adc, intercept)

    println("ADC = ", adc)
    println("intercept = ", intercept)

    return signals, ln_signal, adc, intercept, fit_ln
end

# ============================================================
# Signal processing
# ============================================================

"""
    normalize_by_first_value(signals)

Normalize a vector by its first value, usually S(0).
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

"""
    fit_adc_from_ln_signal(bvals, ln_signal)

Fit the linear model:

    ln(S/S0) = intercept - ADC * b

using a simple least-squares line fit.

Returns:
- `adc`
- `intercept`
"""
function fit_adc_from_ln_signal(bvals, ln_signal)
    valid = .!isnan.(bvals) .& .!isnan.(ln_signal)

    b = bvals[valid]
    y = ln_signal[valid]

    if length(b) < 2
        return NaN, NaN
    end

    b_mean = mean(b)
    y_mean = mean(y)

    slope = sum((b .- b_mean) .* (y .- y_mean)) / sum((b .- b_mean).^2)
    intercept = y_mean - slope * b_mean

    adc = -slope

    return adc, intercept
end

"""
    fitted_ln_signal(bvals, adc, intercept)

Return fitted values for:

    ln(S/S0) = intercept - ADC * b
"""
function fitted_ln_signal(bvals, adc, intercept)
    return intercept .- adc .* bvals
end

# ============================================================
# Result rows and CSV writing
# ============================================================

"""
    signal_header()

Return the CSV header for detailed signal outputs.
"""
function signal_header()
    return [
        "substrate_name",
        "substrate_format",
        "TE",
        "direction_name",
        "gradient_x",
        "gradient_y",
        "gradient_z",
        "b",
        "repeat",
        "nspins",
        "n_inside_spins",
        "n_outside_spins",
        "signal",
        "normalized_signal",
        "ln_S_over_S0",
        "fitted_ln_S_over_S0",
    ]
end

"""
    summary_header()

Return the CSV header for ADC summary outputs.
"""
function summary_header()
    return [
        "substrate_name",
        "substrate_format",
        "TE",
        "direction_name",
        "gradient_x",
        "gradient_y",
        "gradient_z",
        "repeat",
        "nspins",
        "n_inside_spins",
        "n_outside_spins",
        "ADC",
        "intercept",
    ]
end

"""
    append_signal_rows!(rows, ...)

Append detailed per-b-value signal rows to `rows`.
"""
function append_signal_rows!(
    rows,
    substrate_name,
    substrate_format,
    TE,
    direction_name,
    direction,
    bvals,
    repeat_index,
    nspins,
    n_inside,
    n_outside,
    signals,
    ln_signal,
    fit_ln,
)
    normalized = normalize_by_first_value(signals)

    for i in eachindex(bvals)
        push!(
            rows,
            [
                substrate_name,
                substrate_format,
                TE,
                direction_name,
                direction[1],
                direction[2],
                direction[3],
                bvals[i],
                repeat_index,
                nspins,
                n_inside,
                n_outside,
                signals[i],
                normalized[i],
                ln_signal[i],
                fit_ln[i],
            ],
        )
    end
end

"""
    append_summary_row!(rows, ...)

Append one ADC summary row to `rows`.
"""
function append_summary_row!(
    rows,
    substrate_name,
    substrate_format,
    TE,
    direction_name,
    direction,
    repeat_index,
    nspins,
    n_inside,
    n_outside,
    adc,
    intercept,
)
    push!(
        rows,
        [
            substrate_name,
            substrate_format,
            TE,
            direction_name,
            direction[1],
            direction[2],
            direction[3],
            repeat_index,
            nspins,
            n_inside,
            n_outside,
            adc,
            intercept,
        ],
    )
end

"""
    write_rows_csv(output_file, rows)

Write rows to a CSV file using `DelimitedFiles.writedlm`.
"""
function write_rows_csv(output_file::String, rows)
    open(output_file, "w") do io
        writedlm(io, rows, ',')
    end

    println("Saved CSV:")
    println("  ", output_file)
end

"""
    make_output_paths(prefix::String)

Create standard output paths for one simulation run.

Returns:
- `signal_csv`
- `summary_csv`
"""
function make_output_paths(prefix::String)
    output_dir, _ = ensure_output_dirs()

    signal_csv = joinpath(
        output_dir,
        prefix * "_signal.csv",
    )

    summary_csv = joinpath(
        output_dir,
        prefix * "_ADC_summary.csv",
    )

    return signal_csv, summary_csv
end
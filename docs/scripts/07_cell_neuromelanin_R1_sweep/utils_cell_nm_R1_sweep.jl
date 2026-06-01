# ============================================================
# utils_cell_nm_R1_sweep.jl
#
# Utility functions for cell + neuromelanin R1 sweep simulations.
#
# This file is included by scripts in:
#
#     docs/scripts/07_cell_neuromelanin_R1_sweep/
#
# Goal:
# Build a reusable workflow for:
#   - loading a cell SWC substrate;
#   - loading a neuromelanin SWC substrate;
#   - assigning local R1 to the neuromelanin compartment;
#   - keeping global R1/R2 fixed;
#   - running gradient-echo simulations;
#   - extracting compartment-wise signals;
#   - computing NM-MRI contrast ratio:
#
#         CR = (S_i - S_0) / S_0
#
# where S_0 is the baseline signal when R1_NM equals the background R1.
# ============================================================

using MCMRSimulator
using MRIBuilder
using DelimitedFiles
using Random
using Statistics

# ============================================================
# Unit conversion
# ============================================================

"""
    per_second_to_per_ms(rate_s)

Convert a rate from s^-1 to ms^-1.

MCMR simulations use sequence timings in ms, so relaxation rates should
be passed in ms^-1 when used together with TE/TR in ms.

Example:
    R1 = 1 s^-1 = 0.001 ms^-1
"""
function per_second_to_per_ms(rate_s)
    return rate_s / 1000.0
end

"""
    per_ms_to_per_second(rate_ms)

Convert a rate from ms^-1 to s^-1.
"""
function per_ms_to_per_second(rate_ms)
    return rate_ms * 1000.0
end

# ============================================================
# Paths and file handling
# ============================================================

"""
    default_substrate_dir()

Return the folder where example substrate files are stored.
"""
function default_substrate_dir()
    return joinpath(
        "substrate_generation",
        "example_outputs",
    )
end

"""
    ensure_output_dirs()

Create and return the standard output directories.

Returns:
- `output_dir`: docs/data/processed
- `figure_dir`: docs/figures
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
"""
function add_extension_if_missing(filename::String, extension::String)
    if endswith(filename, extension)
        return filename
    else
        return filename * extension
    end
end

"""
    swc_path_from_name(swc_name::String)

Return the relative path to an SWC file stored in the default substrate
directory.
"""
function swc_path_from_name(swc_name::String)
    filename = add_extension_if_missing(swc_name, ".swc")

    return joinpath(
        default_substrate_dir(),
        filename,
    )
end

"""
    check_file_exists(path::String; description = "file")

Check that a file exists and throw a clear error otherwise.
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
    print_selected_files(cell_swc_name, nm_swc_name)

Print the expected cell and neuromelanin SWC files.
"""
function print_selected_files(cell_swc_name::String, nm_swc_name::String)
    println("Selected substrates")
    println("-------------------")
    println("Cell SWC:")
    println("  ", swc_path_from_name(cell_swc_name))
    println("Neuromelanin SWC:")
    println("  ", swc_path_from_name(nm_swc_name))
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
    load_swc_spheres(swc_name::String)

Load sphere positions and radii from an SWC file stored in the default
substrate directory.

The `.swc` extension is optional.
"""
function load_swc_spheres(swc_name::String)
    swc_path = swc_path_from_name(swc_name)

    println("Loading SWC:")
    println("  ", swc_path)

    positions, radii = read_swc_spheres(swc_path)

    println("Loaded successfully.")
    println("Number of spheres: ", length(radii))
    println("Minimum radius: ", minimum(radii))
    println("Maximum radius: ", maximum(radii))

    return positions, radii, swc_path
end

# ============================================================
# Geometry creation
# ============================================================

"""
    create_cell_geometry(positions, radii; local_R1_additional_ms = 0.0)

Create the cell geometry as overlapping spheres.

The cell is the main substrate. In this project, the global R1 is usually
applied to the whole simulation, so the cell-specific local R1 contribution
is often set to zero.

`local_R1_additional_ms` is an additional local R1 contribution in ms^-1.
"""
function create_cell_geometry(
    positions,
    radii;
    local_R1_additional_ms = 0.0,
)
    geometry = Spheres(
        radius = radii,
        position = positions,
        overlapping = true,
        R1_inside = local_R1_additional_ms,
        permeability = 0.0,
    )

    return geometry
end

"""
    create_nm_geometry(positions, radii; local_R1_additional_ms)

Create the neuromelanin geometry as overlapping spheres.

`local_R1_additional_ms` is the additional local R1 contribution in ms^-1.
For example, if the global background R1 is 1 s^-1 and the target total
R1_NM is 50 s^-1, then:

    additional R1_NM = 50 - 1 = 49 s^-1

which is converted to ms^-1 before being passed here.
"""
function create_nm_geometry(
    positions,
    radii;
    local_R1_additional_ms,
)
    geometry = Spheres(
        radius = radii,
        position = positions,
        overlapping = true,
        R1_inside = local_R1_additional_ms,
        permeability = Inf,
    )

    return geometry
end

"""
    create_cell_nm_geometry(cell_positions, cell_radii, nm_positions, nm_radii;
                            R1_background_s, R1_NM_total_s)

Create a two-substrate geometry with:
- substrate 1: cell;
- substrate 2: neuromelanin.

The global/background R1 is applied to the whole simulation.
The neuromelanin geometry receives an additional local R1 contribution so
that its total R1 equals `R1_NM_total_s`.

Returns:
- `geometry`: vector containing cell and NM geometries;
- `R1_NM_additional_s`;
- `R1_NM_additional_ms`.
"""
function create_cell_nm_geometry(
    cell_positions,
    cell_radii,
    nm_positions,
    nm_radii;
    R1_background_s,
    R1_NM_total_s,
)
    R1_NM_additional_s = R1_NM_total_s - R1_background_s

    if R1_NM_additional_s < 0
        @warn """
        R1_NM_total_s is smaller than R1_background_s.

        This gives a negative additional local R1 contribution:
            R1_NM_total_s      = $R1_NM_total_s s^-1
            R1_background_s    = $R1_background_s s^-1
            R1_NM_additional_s = $R1_NM_additional_s s^-1
        """
    end

    R1_NM_additional_ms = per_second_to_per_ms(R1_NM_additional_s)

    cell_geometry = create_cell_geometry(
        cell_positions,
        cell_radii;
        local_R1_additional_ms = 0.0,
    )

    nm_geometry = create_nm_geometry(
        nm_positions,
        nm_radii;
        local_R1_additional_ms = R1_NM_additional_ms,
    )

    geometry = [
        cell_geometry,
        nm_geometry,
    ]

    return geometry, R1_NM_additional_s, R1_NM_additional_ms
end

# ============================================================
# Geometry extent and spin generation
# ============================================================

"""
    sphere_ranges(positions, radii)

Return spatial ranges of a sphere substrate, including radii.

Returns:
- xmin, xmax, ymin, ymax, zmin, zmax
"""
function sphere_ranges(positions, radii)
    xs = [p[1] for p in positions]
    ys = [p[2] for p in positions]
    zs = [p[3] for p in positions]

    xmin = minimum(xs .- radii)
    xmax = maximum(xs .+ radii)

    ymin = minimum(ys .- radii)
    ymax = maximum(ys .+ radii)

    zmin = minimum(zs .- radii)
    zmax = maximum(zs .+ radii)

    return xmin, xmax, ymin, ymax, zmin, zmax
end

"""
    print_sphere_ranges(label, positions, radii)

Print spatial ranges of a sphere substrate.
"""
function print_sphere_ranges(label, positions, radii)
    xmin, xmax, ymin, ymax, zmin, zmax = sphere_ranges(
        positions,
        radii,
    )

    println(label)
    println("  x range: ", xmin, " to ", xmax)
    println("  y range: ", ymin, " to ", ymax)
    println("  z range: ", zmin, " to ", zmax)
end

"""
    combined_bbox_from_two_substrates(cell_positions, cell_radii,
                                      nm_positions, nm_radii; margin = 2.0)

Create a bounding box containing both the cell and neuromelanin substrates.
"""
function combined_bbox_from_two_substrates(
    cell_positions,
    cell_radii,
    nm_positions,
    nm_radii;
    margin = 2.0,
)
    c_xmin, c_xmax, c_ymin, c_ymax, c_zmin, c_zmax =
        sphere_ranges(cell_positions, cell_radii)

    n_xmin, n_xmax, n_ymin, n_ymax, n_zmin, n_zmax =
        sphere_ranges(nm_positions, nm_radii)

    lower = [
        min(c_xmin, n_xmin) - margin,
        min(c_ymin, n_ymin) - margin,
        min(c_zmin, n_zmin) - margin,
    ]

    upper = [
        max(c_xmax, n_xmax) + margin,
        max(c_ymax, n_ymax) + margin,
        max(c_zmax, n_zmax) + margin,
    ]

    return BoundingBox(lower, upper)
end

"""
    create_spins_for_cell_nm(cell_positions, cell_radii,
                             nm_positions, nm_radii;
                             nspins = 50_000,
                             margin = 2.0,
                             seed = 1234)

Create random spins in a bounding box containing both the cell and NM
substrates.

Returns:
- `spins`
- `bbox`
"""
function create_spins_for_cell_nm(
    cell_positions,
    cell_radii,
    nm_positions,
    nm_radii;
    nspins = 50_000,
    margin = 2.0,
    seed = 1234,
)
    Random.seed!(seed)

    bbox = combined_bbox_from_two_substrates(
        cell_positions,
        cell_radii,
        nm_positions,
        nm_radii;
        margin = margin,
    )

    println("Combined bounding box:")
    println(bbox)

    spins = Snapshot(nspins, bbox)

    println("Generated spins: ", length(spins))

    return spins, bbox
end

# ============================================================
# Compartment masks and counts
# ============================================================

"""
    count_inside_geometry(geometry, spins)

Count spins inside a geometry.
"""
function count_inside_geometry(geometry, spins)
    inside_mask = isinside(geometry, spins) .> 0
    return count(inside_mask)
end

"""
    print_compartment_counts(cell_geometry, nm_geometry, spins)

Print spin counts for cell and neuromelanin compartments.
"""
function print_compartment_counts(cell_geometry, nm_geometry, spins)
    inside_cell = count_inside_geometry(cell_geometry, spins)
    inside_nm = count_inside_geometry(nm_geometry, spins)

    println("Compartment spin counts")
    println("-----------------------")
    println("Total spins:      ", length(spins))
    println("Inside cell:      ", inside_cell)
    println("Outside cell:     ", length(spins) - inside_cell)
    println("Inside NM:        ", inside_nm)
    println("Outside NM:       ", length(spins) - inside_nm)

    return inside_cell, inside_nm
end

# ============================================================
# MRI sequence and simulation
# ============================================================

"""
    create_gradient_echo_sequence(; TE = 3, TR = 25)

Create a gradient-echo sequence.

TE and TR are in ms.
"""
function create_gradient_echo_sequence(; TE = 3, TR = 25)
    return GradientEcho(
        TE = TE,
        TR = TR,
        scanner = MRIBuilder.Siemens_Prisma,
    )
end

"""
    create_gradient_echo_simulation(geometry;
                                    TE,
                                    TR,
                                    R1_background_s,
                                    R2_global_s,
                                    diffusivity,
                                    timestep)

Create a GradientEcho simulation with:
- fixed global/background R1;
- fixed global R2;
- local R1 contributions already encoded inside the geometry.

All input rates are given in s^-1 and converted to ms^-1 here.
"""
function create_gradient_echo_simulation(
    geometry;
    TE = 3,
    TR = 25,
    R1_background_s = 1.0,
    R2_global_s = 1.0 / 3.0,
    diffusivity = 1.0,
    timestep = nothing,
)
    sequence = create_gradient_echo_sequence(
        TE = TE,
        TR = TR,
    )

    R1_background_ms = per_second_to_per_ms(R1_background_s)
    R2_global_ms = per_second_to_per_ms(R2_global_s)

    if isnothing(timestep)
        simulation = Simulation(
            sequence,
            diffusivity = diffusivity,
            geometry = geometry,
            R1 = R1_background_ms,
            R2 = R2_global_ms,
        )
    else
        simulation = Simulation(
            sequence,
            diffusivity = diffusivity,
            geometry = geometry,
            R1 = R1_background_ms,
            R2 = R2_global_ms,
            timestep = timestep,
        )
    end

    return simulation
end

"""
    signal_scalar(signal)

Extract a scalar transverse signal magnitude from an MCMR signal object.
"""
function signal_scalar(signal)
    return abs(signal.orient.transverse)
end

"""
    readout_cell_nm_compartments(spins, simulation; skip_TR = 20)

Run a GradientEcho readout and extract signals from relevant compartments.

Substrate convention:
- geometry_index = 1: cell
- geometry_index = 2: neuromelanin

Returns a named tuple with:
- total
- inside_cell
- outside_cell
- inside_nm
- outside_nm
"""
function readout_cell_nm_compartments(
    spins,
    simulation;
    skip_TR = 20,
)
    compartment_signals = readout(
        spins,
        simulation,
        skip_TR = skip_TR,
        subset = [
            Subset(),                              # total
            Subset(inside = true, geometry_index = 1),
            Subset(inside = false, geometry_index = 1),
            Subset(inside = true, geometry_index = 2),
            Subset(inside = false, geometry_index = 2),
        ],
    )

    return (
        total = signal_scalar(compartment_signals[1]),
        inside_cell = signal_scalar(compartment_signals[2]),
        outside_cell = signal_scalar(compartment_signals[3]),
        inside_nm = signal_scalar(compartment_signals[4]),
        outside_nm = signal_scalar(compartment_signals[5]),
        n_total = compartment_signals[1].nspins,
        n_inside_cell = compartment_signals[2].nspins,
        n_outside_cell = compartment_signals[3].nspins,
        n_inside_nm = compartment_signals[4].nspins,
        n_outside_nm = compartment_signals[5].nspins,
    )
end

"""
    run_one_R1_NM_simulation(spins, cell_positions, cell_radii,
                             nm_positions, nm_radii; kwargs...)

Create a cell+NM geometry for one R1_NM value, run GradientEcho,
and return compartment-wise signals.
"""
function run_one_R1_NM_simulation(
    spins,
    cell_positions,
    cell_radii,
    nm_positions,
    nm_radii;
    R1_background_s,
    R1_NM_total_s,
    R2_global_s,
    TE,
    TR,
    diffusivity,
    timestep,
    skip_TR,
)
    geometry, R1_NM_additional_s, R1_NM_additional_ms =
        create_cell_nm_geometry(
            cell_positions,
            cell_radii,
            nm_positions,
            nm_radii;
            R1_background_s = R1_background_s,
            R1_NM_total_s = R1_NM_total_s,
        )

    simulation = create_gradient_echo_simulation(
        geometry;
        TE = TE,
        TR = TR,
        R1_background_s = R1_background_s,
        R2_global_s = R2_global_s,
        diffusivity = diffusivity,
        timestep = timestep,
    )

    signals = readout_cell_nm_compartments(
        spins,
        simulation;
        skip_TR = skip_TR,
    )

    return signals, R1_NM_additional_s, R1_NM_additional_ms, geometry
end

# ============================================================
# Contrast ratio
# ============================================================

"""
    contrast_ratio(signal, baseline_signal)

Compute contrast ratio:

    CR = (signal - baseline_signal) / baseline_signal
"""
function contrast_ratio(signal, baseline_signal)
    if isnan(signal) || isnan(baseline_signal) || baseline_signal == 0
        return NaN
    end

    return (signal - baseline_signal) / baseline_signal
end

"""
    compute_CRs(signals, baseline_signals)

Compute contrast ratios for all signal compartments.
"""
function compute_CRs(signals, baseline_signals)
    return (
        CR_total = contrast_ratio(signals.total, baseline_signals.total),
        CR_inside_cell = contrast_ratio(signals.inside_cell, baseline_signals.inside_cell),
        CR_outside_cell = contrast_ratio(signals.outside_cell, baseline_signals.outside_cell),
        CR_inside_nm = contrast_ratio(signals.inside_nm, baseline_signals.inside_nm),
        CR_outside_nm = contrast_ratio(signals.outside_nm, baseline_signals.outside_nm),
    )
end

# ============================================================
# CSV output
# ============================================================

"""
    R1_sweep_header()

Return CSV header for the R1 sweep simulation output.
"""
function R1_sweep_header()
    return [
        "cell_swc_name",
        "nm_swc_name",
        "R1_background_s",
        "R1_background_ms",
        "R1_NM_total_s",
        "R1_NM_total_ms",
        "R1_NM_additional_s",
        "R1_NM_additional_ms",
        "R2_global_s",
        "R2_global_ms",
        "TE",
        "TR",
        "diffusivity",
        "timestep_ms",
        "skip_TR",
        "n_total",
        "n_inside_cell",
        "n_outside_cell",
        "n_inside_nm",
        "n_outside_nm",
        "total_signal",
        "inside_cell_signal",
        "outside_cell_signal",
        "inside_nm_signal",
        "outside_nm_signal",
        "CR_total",
        "CR_inside_cell",
        "CR_outside_cell",
        "CR_inside_nm",
        "CR_outside_nm",
    ]
end

"""
    append_R1_sweep_row!(rows, ...)

Append one row to the R1 sweep output table.
"""
function append_R1_sweep_row!(
    rows,
    cell_swc_name,
    nm_swc_name,
    R1_background_s,
    R1_NM_total_s,
    R1_NM_additional_s,
    R1_NM_additional_ms,
    R2_global_s,
    TE,
    TR,
    diffusivity,
    timestep,
    skip_TR,
    signals,
    CRs,
)
    R1_background_ms = per_second_to_per_ms(R1_background_s)
    R1_NM_total_ms = per_second_to_per_ms(R1_NM_total_s)
    R2_global_ms = per_second_to_per_ms(R2_global_s)

    push!(
        rows,
        [
            cell_swc_name,
            nm_swc_name,
            R1_background_s,
            R1_background_ms,
            R1_NM_total_s,
            R1_NM_total_ms,
            R1_NM_additional_s,
            R1_NM_additional_ms,
            R2_global_s,
            R2_global_ms,
            TE,
            TR,
            diffusivity,
            isnothing(timestep) ? "automatic" : timestep,
            skip_TR,
            signals.n_total,
            signals.n_inside_cell,
            signals.n_outside_cell,
            signals.n_inside_nm,
            signals.n_outside_nm,
            signals.total,
            signals.inside_cell,
            signals.outside_cell,
            signals.inside_nm,
            signals.outside_nm,
            CRs.CR_total,
            CRs.CR_inside_cell,
            CRs.CR_outside_cell,
            CRs.CR_inside_nm,
            CRs.CR_outside_nm,
        ],
    )
end

"""
    local_R1_check_header()

Return CSV header for local R1 check output.
"""
function local_R1_check_header()
    return [
        "cell_swc_name",
        "nm_swc_name",
        "R1_background_s",
        "R1_NM_total_s",
        "R1_NM_additional_s",
        "x",
        "y",
        "z",
        "inside_cell",
        "inside_nm",
        "expected_R1_s",
    ]
end

"""
    append_local_R1_check_rows!(rows, ...)

Append simple local R1 check rows based on sphere center positions.

This is a lightweight diagnostic. It checks whether chosen points are
inside the cell and/or inside the NM geometry and records the expected R1.
"""
function append_local_R1_check_rows!(
    rows,
    cell_swc_name,
    nm_swc_name,
    R1_background_s,
    R1_NM_total_s,
    R1_NM_additional_s,
    cell_geometry,
    nm_geometry,
    check_positions,
)
    for p in check_positions
        snapshot = Snapshot([p])

        inside_cell = (isinside(cell_geometry, snapshot)[1] > 0)
        inside_nm = (isinside(nm_geometry, snapshot)[1] > 0)

        expected_R1_s = inside_nm ? R1_NM_total_s : R1_background_s

        push!(
            rows,
            [
                cell_swc_name,
                nm_swc_name,
                R1_background_s,
                R1_NM_total_s,
                R1_NM_additional_s,
                p[1],
                p[2],
                p[3],
                inside_cell,
                inside_nm,
                expected_R1_s,
            ],
        )
    end
end

"""
    write_rows_csv(output_file, rows)

Write rows to a CSV file.
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

Return standard output paths for the R1 sweep.

Returns:
- `sweep_csv`
- `local_R1_check_csv`
"""
function make_output_paths(prefix::String)
    output_dir, _ = ensure_output_dirs()

    sweep_csv = joinpath(
        output_dir,
        prefix * "_R1_sweep.csv",
    )

    local_R1_check_csv = joinpath(
        output_dir,
        prefix * "_local_R1_check.csv",
    )

    return sweep_csv, local_R1_check_csv
end

# ============================================================
# R1 value grid
# ============================================================

"""
    default_R1_NM_values_s()

Return a default R1_NM sweep in s^-1.

This includes:
- baseline 1 s^-1;
- values from 10 to 100 s^-1 in steps of 10;
- values from 200 to 1000 s^-1 in steps of 100.

This follows the supervisor's suggestion and keeps 50 s^-1 included.
"""
function default_R1_NM_values_s()
    return vcat(
        [1.0],
        collect(10.0:10.0:100.0),
        collect(200.0:100.0:1000.0),
    )
end

"""
    coarse_R1_NM_values_s()

Return a coarser R1_NM sweep in s^-1.

Useful if the full sweep takes too long.
"""
function coarse_R1_NM_values_s()
    return vcat(
        [1.0],
        collect(20.0:20.0:100.0),
        collect(200.0:200.0:1000.0),
    )
end
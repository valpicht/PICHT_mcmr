# ============================================================
# utils_cell_neuromelanin.jl
#
# Utility functions for two-substrate simulations:
#
#   1. cell substrate;
#   2. neuromelanin substrate.
#
# Both substrates are loaded from SWC files and converted into
# overlapping-sphere geometries.
#
# This file is included by scripts in:
#
#     docs/scripts/06_two_substrate_cell_neuromelanin/
#
# The goal is to keep the main scripts short and reusable.
# ============================================================

using MCMRSimulator
using MRIBuilder
using DelimitedFiles
using Random
using Statistics

# ============================================================
# Paths and file handling
# ============================================================

"""
    default_substrate_dir()

Return the folder where example substrate files are stored.

Expected location:

    substrate_generation/example_outputs/
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
    add_swc_extension_if_missing(filename::String)

Add `.swc` to a file name if it is missing.
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

Return the full relative path to an SWC file stored in the default
substrate directory.
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
    print_selected_substrates(cell_swc_name, nm_swc_name)

Print the selected cell and neuromelanin substrate paths.
"""
function print_selected_substrates(cell_swc_name::String, nm_swc_name::String)
    println("Selected substrates")
    println("-------------------")

    println("Cell SWC name:")
    println("  ", cell_swc_name)

    println("Cell SWC path:")
    println("  ", swc_path_from_name(cell_swc_name))

    println("Neuromelanin SWC name:")
    println("  ", nm_swc_name)

    println("Neuromelanin SWC path:")
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
- `positions`: vector of 3-element position vectors `[x, y, z]`
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
    load_swc_spheres(swc_name::String)

Load an SWC file from `substrate_generation/example_outputs/`.

Returns:
- `positions`
- `radii`
- `swc_path`
"""
function load_swc_spheres(swc_name::String)
    swc_path = swc_path_from_name(swc_name)

    println("Loading SWC:")
    println("  ", swc_path)

    positions, radii = read_swc_spheres(swc_path)

    println("SWC loaded successfully.")
    println("Number of spheres: ", length(radii))
    println("Minimum radius: ", minimum(radii))
    println("Maximum radius: ", maximum(radii))

    return positions, radii, swc_path
end

# ============================================================
# Geometry creation
# ============================================================

"""
    create_cell_geometry(positions, radii; cell_R1_inside = 0.0, permeability = 0.0)

Create an overlapping-sphere geometry for the cell substrate.

`cell_R1_inside` is an additional local R1 contribution assigned inside
the cell geometry.

`permeability = 0.0` keeps the cell boundary impermeable in the first
validation scripts.
"""
function create_cell_geometry(
    positions,
    radii;
    cell_R1_inside = 0.0,
    permeability = 0.0,
)
    geometry = Spheres(
        radius = radii,
        position = positions,
        overlapping = true,
        R1_inside = cell_R1_inside,
        permeability = permeability,
    )

    println("Cell geometry created.")
    println("Geometry type: ", typeof(geometry))
    println("Cell R1_inside: ", cell_R1_inside)
    println("Cell permeability: ", permeability)

    return geometry
end

"""
    create_neuromelanin_geometry(positions, radii; nm_R1_inside = 0.0, permeability = Inf)

Create an overlapping-sphere geometry for the neuromelanin substrate.

`nm_R1_inside` is an additional local R1 contribution assigned inside
the neuromelanin geometry.

`permeability = Inf` was used in the project to allow spins to move
through the neuromelanin region while still assigning local relaxation
properties.
"""
function create_neuromelanin_geometry(
    positions,
    radii;
    nm_R1_inside = 0.0,
    permeability = Inf,
)
    geometry = Spheres(
        radius = radii,
        position = positions,
        overlapping = true,
        R1_inside = nm_R1_inside,
        permeability = permeability,
    )

    println("Neuromelanin geometry created.")
    println("Geometry type: ", typeof(geometry))
    println("Neuromelanin R1_inside: ", nm_R1_inside)
    println("Neuromelanin permeability: ", permeability)

    return geometry
end

"""
    load_cell_and_neuromelanin_geometries(cell_swc_name, nm_swc_name; kwargs...)

Load the cell and neuromelanin SWC files and convert both into
overlapping-sphere geometries.

Returns:
- `geometry`: vector containing `[cell_geometry, nm_geometry]`
- `cell_positions`
- `cell_radii`
- `nm_positions`
- `nm_radii`
"""
function load_cell_and_neuromelanin_geometries(
    cell_swc_name::String,
    nm_swc_name::String;
    cell_R1_inside = 0.0,
    nm_R1_inside = 0.0,
    cell_permeability = 0.0,
    nm_permeability = Inf,
)
    println("Loading cell substrate...")
    cell_positions, cell_radii, _ = load_swc_spheres(cell_swc_name)

    println("\nLoading neuromelanin substrate...")
    nm_positions, nm_radii, _ = load_swc_spheres(nm_swc_name)

    println("\nCreating cell geometry...")
    cell_geometry = create_cell_geometry(
        cell_positions,
        cell_radii;
        cell_R1_inside = cell_R1_inside,
        permeability = cell_permeability,
    )

    println("\nCreating neuromelanin geometry...")
    nm_geometry = create_neuromelanin_geometry(
        nm_positions,
        nm_radii;
        nm_R1_inside = nm_R1_inside,
        permeability = nm_permeability,
    )

    geometry = [
        cell_geometry,
        nm_geometry,
    ]

    return geometry, cell_positions, cell_radii, nm_positions, nm_radii
end

# ============================================================
# Geometry inspection and bounding boxes
# ============================================================

"""
    sphere_ranges(positions, radii)

Return spatial ranges of a set of spheres, including their radii.

Returns:
- `xmin, xmax, ymin, ymax, zmin, zmax`
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

Print the spatial extent of a sphere-based substrate.
"""
function print_sphere_ranges(label::String, positions, radii)
    xmin, xmax, ymin, ymax, zmin, zmax = sphere_ranges(positions, radii)

    println(label)
    println("x range: ", xmin, " to ", xmax)
    println("y range: ", ymin, " to ", ymax)
    println("z range: ", zmin, " to ", zmax)
end

"""
    make_bbox_from_two_substrates(cell_positions, cell_radii, nm_positions, nm_radii; margin = 2.0)

Create a bounding box that contains both the cell and neuromelanin
substrates.
"""
function make_bbox_from_two_substrates(
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

    xmin = min(c_xmin, n_xmin)
    xmax = max(c_xmax, n_xmax)

    ymin = min(c_ymin, n_ymin)
    ymax = max(c_ymax, n_ymax)

    zmin = min(c_zmin, n_zmin)
    zmax = max(c_zmax, n_zmax)

    return BoundingBox(
        [xmin - margin, ymin - margin, zmin - margin],
        [xmax + margin, ymax + margin, zmax + margin],
    )
end

"""
    create_spins_for_two_substrates(...; nspins = 10_000, margin = 2.0, seed = 1234)

Create random spins in a bounding box containing both substrates.

Returns:
- `spins`
- `bbox`
"""
function create_spins_for_two_substrates(
    cell_positions,
    cell_radii,
    nm_positions,
    nm_radii;
    nspins = 10_000,
    margin = 2.0,
    seed = 1234,
)
    Random.seed!(seed)

    bbox = make_bbox_from_two_substrates(
        cell_positions,
        cell_radii,
        nm_positions,
        nm_radii;
        margin = margin,
    )

    println("Bounding box:")
    println(bbox)

    spins = Snapshot(nspins, bbox)

    println("Number of generated spins: ", length(spins))

    return spins, bbox
end

# ============================================================
# Inside/outside checks
# ============================================================

"""
    inside_mask(geometry, spins)

Return a boolean mask indicating which spins are inside a geometry.
"""
function inside_mask(geometry, spins)
    return isinside(geometry, spins) .> 0
end

"""
    count_inside_outside(label, geometry, spins)

Print and return the number of spins inside and outside a geometry.
"""
function count_inside_outside(label::String, geometry, spins)
    mask = inside_mask(geometry, spins)

    n_inside = count(mask)
    n_outside = length(spins) - n_inside

    println(label)
    println("  Inside spins: ", n_inside)
    println("  Outside spins: ", n_outside)

    return n_inside, n_outside
end

"""
    count_cell_nm_compartments(geometry, spins)

Count spins inside/outside the cell and inside/outside neuromelanin.

Assumes:
- `geometry[1]` is the cell geometry;
- `geometry[2]` is the neuromelanin geometry.
"""
function count_cell_nm_compartments(geometry, spins)
    cell_geometry = geometry[1]
    nm_geometry = geometry[2]

    println("Compartment spin counts")
    println("-----------------------")

    n_cell_inside, n_cell_outside = count_inside_outside(
        "Cell geometry:",
        cell_geometry,
        spins,
    )

    n_nm_inside, n_nm_outside = count_inside_outside(
        "Neuromelanin geometry:",
        nm_geometry,
        spins,
    )

    return n_cell_inside, n_cell_outside, n_nm_inside, n_nm_outside
end

# ============================================================
# Expected local R1 checks
# ============================================================

"""
    expected_local_R1_at_position(position, geometry; global_R1 = 0.0, cell_R1_inside = 0.0, nm_R1_inside = 0.0)

Compute the expected local R1 at one position based on whether the
position lies inside the cell and/or neuromelanin geometries.

This function is a practical project-level check of the intended R1
assignment.

Assumes:
- `geometry[1]` is the cell geometry;
- `geometry[2]` is the neuromelanin geometry.

The expected value is:

    global_R1
    + cell_R1_inside if inside cell
    + nm_R1_inside if inside neuromelanin
"""
function expected_local_R1_at_position(
    position,
    geometry;
    global_R1 = 0.0,
    cell_R1_inside = 0.0,
    nm_R1_inside = 0.0,
)
    single_spin = Snapshot([position])

    inside_cell = isinside(geometry[1], single_spin)[1] > 0
    inside_nm = isinside(geometry[2], single_spin)[1] > 0

    expected_R1 = global_R1

    if inside_cell
        expected_R1 += cell_R1_inside
    end

    if inside_nm
        expected_R1 += nm_R1_inside
    end

    return expected_R1, inside_cell, inside_nm
end

"""
    build_R1_check_rows(test_points, geometry; kwargs...)

Build CSV rows for checking expected local R1 at selected positions.

`test_points` should be a vector of pairs:

    [
        ("label", [x, y, z]),
        ...
    ]
"""
function build_R1_check_rows(
    test_points,
    geometry;
    global_R1 = 0.0,
    cell_R1_inside = 0.0,
    nm_R1_inside = 0.0,
)
    rows = Any[]

    push!(
        rows,
        [
            "label",
            "x",
            "y",
            "z",
            "inside_cell",
            "inside_neuromelanin",
            "global_R1",
            "cell_R1_inside",
            "nm_R1_inside",
            "expected_local_R1",
        ],
    )

    for (label, position) in test_points
        expected_R1, inside_cell, inside_nm =
            expected_local_R1_at_position(
                position,
                geometry;
                global_R1 = global_R1,
                cell_R1_inside = cell_R1_inside,
                nm_R1_inside = nm_R1_inside,
            )

        push!(
            rows,
            [
                label,
                position[1],
                position[2],
                position[3],
                inside_cell,
                inside_nm,
                global_R1,
                cell_R1_inside,
                nm_R1_inside,
                expected_R1,
            ],
        )
    end

    return rows
end

"""
    default_R1_test_points(cell_positions, nm_positions)

Create simple default test points for local R1 checks.

Uses:
- first cell sphere center;
- first neuromelanin sphere center;
- a point far outside both geometries.
"""
function default_R1_test_points(cell_positions, nm_positions)
    cell_center = cell_positions[1]
    nm_center = nm_positions[1]

    outside_point = [
        maximum([cell_center[1], nm_center[1]]) + 100.0,
        maximum([cell_center[2], nm_center[2]]) + 100.0,
        maximum([cell_center[3], nm_center[3]]) + 100.0,
    ]

    return [
        ("cell_first_sphere_center", cell_center),
        ("neuromelanin_first_sphere_center", nm_center),
        ("outside_point", outside_point),
    ]
end

# ============================================================
# MRI sequence and simulation
# ============================================================

"""
    create_gradient_echo_sequence(; TE = 3, TR = 25)

Create a gradient-echo sequence using MRIBuilder.
"""
function create_gradient_echo_sequence(; TE = 3, TR = 25)
    return GradientEcho(
        TE = TE,
        TR = TR,
        scanner = MRIBuilder.Siemens_Prisma,
    )
end

"""
    create_gradient_echo_simulation(geometry; TE, TR, diffusivity, global_R1, global_R2)

Create a gradient-echo simulation for the combined cell + neuromelanin
geometry.

`geometry` should usually be:

    [cell_geometry, nm_geometry]
"""
function create_gradient_echo_simulation(
    geometry;
    TE = 3,
    TR = 25,
    diffusivity = 1.0,
    global_R1 = 0.0,
    global_R2 = 0.0,
)
    sequence = create_gradient_echo_sequence(
        TE = TE,
        TR = TR,
    )

    simulation = Simulation(
        sequence,
        diffusivity = diffusivity,
        geometry = geometry,
        R1 = global_R1,
        R2 = global_R2,
    )

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
    run_gradient_echo_compartment_readout(spins, simulation; skip_TR = 20)

Run a gradient-echo readout and extract compartment-wise signals.

Returns a named tuple with:
- total
- inside_cell
- outside_cell
- inside_neuromelanin

The subsets assume:
- geometry_index = 1 is the cell;
- geometry_index = 2 is neuromelanin.
"""
function run_gradient_echo_compartment_readout(
    spins,
    simulation;
    skip_TR = 20,
)
    total_signal = readout(
        spins,
        simulation;
        skip_TR = skip_TR,
    )

    compartment_signals = readout(
        spins,
        simulation;
        skip_TR = skip_TR,
        subset = [
            Subset(inside = true, geometry_index = 1),
            Subset(inside = false, geometry_index = 1),
            Subset(inside = true, geometry_index = 2),
        ],
    )

    inside_cell_signal = compartment_signals[1]
    outside_cell_signal = compartment_signals[2]
    inside_nm_signal = compartment_signals[3]

    return (
        total = signal_scalar(total_signal),
        inside_cell = signal_scalar(inside_cell_signal),
        outside_cell = signal_scalar(outside_cell_signal),
        inside_neuromelanin = signal_scalar(inside_nm_signal),
        total_nspins = total_signal.nspins,
        inside_cell_nspins = inside_cell_signal.nspins,
        outside_cell_nspins = outside_cell_signal.nspins,
        inside_neuromelanin_nspins = inside_nm_signal.nspins,
    )
end

# ============================================================
# CSV helpers
# ============================================================

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
    R1_check_output_path(prefix)

Return a standard output path for local R1 check results.
"""
function R1_check_output_path(prefix::String)
    output_dir, _ = ensure_output_dirs()

    return joinpath(
        output_dir,
        prefix * "_local_R1_check.csv",
    )
end

"""
    gradient_echo_signal_output_path(prefix)

Return a standard output path for gradient-echo compartment signal results.
"""
function gradient_echo_signal_output_path(prefix::String)
    output_dir, _ = ensure_output_dirs()

    return joinpath(
        output_dir,
        prefix * "_gradient_echo_compartment_signals.csv",
    )
end

"""
    gradient_echo_signal_header()

Return the CSV header for gradient-echo compartment signal outputs.
"""
function gradient_echo_signal_header()
    return [
        "cell_swc_name",
        "neuromelanin_swc_name",
        "TE",
        "TR",
        "diffusivity",
        "global_R1",
        "global_R2",
        "cell_R1_inside",
        "nm_R1_inside",
        "nspins",
        "skip_TR",
        "total_signal",
        "inside_cell_signal",
        "outside_cell_signal",
        "inside_neuromelanin_signal",
        "total_nspins",
        "inside_cell_nspins",
        "outside_cell_nspins",
        "inside_neuromelanin_nspins",
    ]
end

"""
    append_gradient_echo_signal_row!(rows, ...)

Append one row of gradient-echo compartment signal results.
"""
function append_gradient_echo_signal_row!(
    rows,
    cell_swc_name,
    nm_swc_name,
    TE,
    TR,
    diffusivity,
    global_R1,
    global_R2,
    cell_R1_inside,
    nm_R1_inside,
    nspins,
    skip_TR,
    result,
)
    push!(
        rows,
        [
            cell_swc_name,
            nm_swc_name,
            TE,
            TR,
            diffusivity,
            global_R1,
            global_R2,
            cell_R1_inside,
            nm_R1_inside,
            nspins,
            skip_TR,
            result.total,
            result.inside_cell,
            result.outside_cell,
            result.inside_neuromelanin,
            result.total_nspins,
            result.inside_cell_nspins,
            result.outside_cell_nspins,
            result.inside_neuromelanin_nspins,
        ],
    )
end
# ============================================================
# 01_check_two_substrate_geometry_and_R1.jl
#
# Goal:
# Check the two-substrate cell + neuromelanin setup before running
# the full gradient-echo simulation.
#
# This script:
#   1. loads the cell SWC substrate;
#   2. loads the neuromelanin SWC substrate;
#   3. creates two overlapping-sphere geometries;
#   4. assigns local R1 contributions to the cell and NM regions;
#   5. creates a common bounding box and random spins;
#   6. counts inside/outside spins for each geometry;
#   7. checks the expected local R1 at selected positions;
#   8. saves a local R1 check CSV.
#
# Run from the repository root with:
#
#   julia +1.11 --project=. docs/scripts/06_two_substrate_cell_neuromelanin/01_check_two_substrate_geometry_and_R1.jl
# ============================================================

using MCMRSimulator
using MRIBuilder
using DelimitedFiles
using Random
using Statistics

include(joinpath(@__DIR__, "utils_cell_neuromelanin.jl"))

# ============================================================
# User parameters
# ============================================================

# Input SWC names, without the .swc extension.
#
# Files must be located in:
#
#   substrate_generation/example_outputs/
#
const CELL_SWC_NAME = "cell_soma_r19_one_branche_neurons_list_new"
const NM_SWC_NAME = "neuromelanin_sphere_r15_neurons_list_new"

# Output prefix.
#
# The script will create:
#
#   docs/data/processed/<OUTPUT_PREFIX>_local_R1_check.csv
#
const OUTPUT_PREFIX = "cell_r19_nm_r15"

# Spin generation parameters.
const NSPINS = 10_000
const BBOX_MARGIN = 2.0
const RANDOM_SEED = 1234

# Local relaxation parameters.
#
# IMPORTANT:
# MCMR uses time in ms in the simulation scripts.
# Therefore R1 values here should be in ms^-1.
#
# Example:
#   T1 = 1000 ms  -> R1 = 1 / 1000 ms^-1 = 0.001 ms^-1
#
# In the project, these local R1 values are used as additional local
# R1 contributions inside the corresponding geometry.
#
const GLOBAL_R1 = 0.0

# For the check script, use clearly different values so that the local
# R1 assignment is easy to verify.
const CELL_R1_INSIDE = 0.0005      # ms^-1
const NM_R1_INSIDE = 0.0010        # ms^-1

# Permeability settings used in the project tests.
#
# The cell is kept impermeable in first validation scripts.
# The neuromelanin region is assigned local relaxation properties while
# allowing spins to pass through it.
#
const CELL_PERMEABILITY = 0.0
const NM_PERMEABILITY = Inf

# ============================================================
# Main script
# ============================================================

function main()
    println("============================================================")
    println("01 — Check two-substrate geometry and local R1")
    println("============================================================")

    output_dir, _ = ensure_output_dirs()

    println("\nConfiguration")
    println("-------------")
    println("cell SWC:              ", CELL_SWC_NAME)
    println("neuromelanin SWC:      ", NM_SWC_NAME)
    println("output prefix:         ", OUTPUT_PREFIX)
    println("nspins:                ", NSPINS)
    println("bbox margin:           ", BBOX_MARGIN)
    println("random seed:           ", RANDOM_SEED)
    println("global R1:             ", GLOBAL_R1, " ms^-1")
    println("cell R1_inside:        ", CELL_R1_INSIDE, " ms^-1")
    println("neuromelanin R1_inside:", NM_R1_INSIDE, " ms^-1")
    println("cell permeability:     ", CELL_PERMEABILITY)
    println("NM permeability:       ", NM_PERMEABILITY)

    println("\nSelected input files")
    println("--------------------")
    print_selected_substrates(
        CELL_SWC_NAME,
        NM_SWC_NAME,
    )

    println("\nLoading geometries...")
    geometry, cell_positions, cell_radii, nm_positions, nm_radii =
        load_cell_and_neuromelanin_geometries(
            CELL_SWC_NAME,
            NM_SWC_NAME;
            cell_R1_inside = CELL_R1_INSIDE,
            nm_R1_inside = NM_R1_INSIDE,
            cell_permeability = CELL_PERMEABILITY,
            nm_permeability = NM_PERMEABILITY,
        )

    println("\nGeometry order")
    println("--------------")
    println("geometry[1] = cell")
    println("geometry[2] = neuromelanin")

    println("\nSubstrate extents")
    println("-----------------")
    print_sphere_ranges(
        "Cell substrate:",
        cell_positions,
        cell_radii,
    )

    print_sphere_ranges(
        "Neuromelanin substrate:",
        nm_positions,
        nm_radii,
    )

    println("\nCreating spins in common bounding box...")
    spins, bbox = create_spins_for_two_substrates(
        cell_positions,
        cell_radii,
        nm_positions,
        nm_radii;
        nspins = NSPINS,
        margin = BBOX_MARGIN,
        seed = RANDOM_SEED,
    )

    println("\nCompartment spin counts")
    println("-----------------------")
    count_cell_nm_compartments(
        geometry,
        spins,
    )

    println("\nBuilding local R1 check points...")
    test_points = default_R1_test_points(
        cell_positions,
        nm_positions,
    )

    println("Test points:")
    for (label, position) in test_points
        println("  ", label, " -> ", position)
    end

    println("\nComputing expected local R1 values...")
    R1_rows = build_R1_check_rows(
        test_points,
        geometry;
        global_R1 = GLOBAL_R1,
        cell_R1_inside = CELL_R1_INSIDE,
        nm_R1_inside = NM_R1_INSIDE,
    )

    R1_output_csv = R1_check_output_path(OUTPUT_PREFIX)

    write_rows_csv(
        R1_output_csv,
        R1_rows,
    )

    println("\nLocal R1 check results")
    println("----------------------")
    for row in R1_rows[2:end]
        label = row[1]
        inside_cell = row[5]
        inside_nm = row[6]
        expected_R1 = row[10]

        println("Point: ", label)
        println("  inside cell:        ", inside_cell)
        println("  inside neuromelanin:", inside_nm)
        println("  expected local R1:  ", expected_R1, " ms^-1")
    end

    println("\nInterpretation")
    println("--------------")
    println("The expected local R1 is computed as:")
    println("")
    println("    global_R1")
    println("    + cell_R1_inside if the point is inside the cell geometry")
    println("    + nm_R1_inside if the point is inside the neuromelanin geometry")
    println("")
    println("This script checks the intended local R1 assignment before running")
    println("the full gradient-echo simulation.")

    println("\nCheck completed successfully.")
    println("Output CSV:")
    println("  ", R1_output_csv)
end

main()
# ============================================================
# 02_run_cell_neuromelanin_gradient_echo.jl
#
# Goal:
# Run a gradient-echo simulation with two substrates:
#
#   1. cell substrate;
#   2. neuromelanin substrate.
#
# Both substrates are loaded from SWC files and converted into
# overlapping-sphere geometries.
#
# This script:
#   1. loads the cell SWC substrate;
#   2. loads the neuromelanin SWC substrate;
#   3. creates two overlapping-sphere geometries;
#   4. assigns local R1 contributions;
#   5. sets a global R2;
#   6. generates spins in a common bounding box;
#   7. runs a GradientEcho simulation;
#   8. extracts total and compartment-wise signals;
#   9. saves the results to CSV.
#
# Run from the repository root with:
#
#   julia +1.11 --project=. docs/scripts/06_two_substrate_cell_neuromelanin/02_run_cell_neuromelanin_gradient_echo.jl
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
#   docs/data/processed/<OUTPUT_PREFIX>_gradient_echo_compartment_signals.csv
#
const OUTPUT_PREFIX = "cell_r19_nm_r15"

# Spin generation parameters.
const NSPINS = 50_000
const BBOX_MARGIN = 2.0
const RANDOM_SEED = 1234

# Gradient-echo sequence parameters.
const TE = 3       # ms
const TR = 25      # ms

# Simulation parameters.
const DIFFUSIVITY = 1.0
const SKIP_TR = 20

# Relaxation parameters.
#
# IMPORTANT:
# In these scripts, time is expressed in ms.
# Therefore R1 and R2 values should be in ms^-1.
#
# Examples:
#   T1 = 1000 ms  -> R1 = 1 / 1000 ms^-1 = 0.001 ms^-1
#   T2 = 3 ms     -> R2 = 1 / 3 ms^-1 ≈ 0.333 ms^-1
#
# GLOBAL_R1 is applied everywhere.
# CELL_R1_INSIDE and NM_R1_INSIDE are additional local R1 contributions
# inside the corresponding geometries.
#
const GLOBAL_R1 = 0.0

# In the project, the cell R1 was often kept fixed, while NM R1 was varied.
# These values are additional local R1 contributions in ms^-1.
const CELL_R1_INSIDE = 0.0005

# Example NM R1 values to compare.
#
# You can change this list to test different NM relaxation conditions.
# For example:
#   0.0005 ms^-1 corresponds to T1 = 2000 ms
#   0.0010 ms^-1 corresponds to T1 = 1000 ms
#   1.0000 ms^-1 corresponds to T1 = 1 ms
#
const NM_R1_INSIDE_VALUES = [
    0.0005,
    0.0010,
    1.0000,
]

# Global R2.
#
# Example used during the project:
#   R2 = 1 / 3 ms^-1
#
const GLOBAL_R2 = 1.0 / 3.0

# Permeability settings.
#
# The cell boundary is kept impermeable in the first validation scripts.
# Neuromelanin is assigned local relaxation while allowing spins to pass
# through it.
#
const CELL_PERMEABILITY = 0.0
const NM_PERMEABILITY = Inf

# ============================================================
# Helper
# ============================================================

"""
    run_one_nm_R1_condition(nm_R1_inside)

Run the full two-substrate gradient-echo simulation for one
neuromelanin R1 value.

Returns:
- one named tuple containing compartment signals and spin counts.
"""
function run_one_nm_R1_condition(nm_R1_inside)
    println("\n============================================================")
    println("Running NM R1 condition")
    println("============================================================")
    println("NM R1_inside = ", nm_R1_inside, " ms^-1")

    println("\nLoading geometries...")
    geometry, cell_positions, cell_radii, nm_positions, nm_radii =
        load_cell_and_neuromelanin_geometries(
            CELL_SWC_NAME,
            NM_SWC_NAME;
            cell_R1_inside = CELL_R1_INSIDE,
            nm_R1_inside = nm_R1_inside,
            cell_permeability = CELL_PERMEABILITY,
            nm_permeability = NM_PERMEABILITY,
        )

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

    println("\nCreating GradientEcho simulation...")
    simulation = create_gradient_echo_simulation(
        geometry;
        TE = TE,
        TR = TR,
        diffusivity = DIFFUSIVITY,
        global_R1 = GLOBAL_R1,
        global_R2 = GLOBAL_R2,
    )

    println("\nRunning readout...")
    result = @time run_gradient_echo_compartment_readout(
        spins,
        simulation;
        skip_TR = SKIP_TR,
    )

    println("\nSignals")
    println("-------")
    println("total signal:                 ", result.total)
    println("inside cell signal:           ", result.inside_cell)
    println("outside cell signal:          ", result.outside_cell)
    println("inside neuromelanin signal:   ", result.inside_neuromelanin)

    println("\nSpin counts from readout")
    println("------------------------")
    println("total nspins:                 ", result.total_nspins)
    println("inside cell nspins:           ", result.inside_cell_nspins)
    println("outside cell nspins:          ", result.outside_cell_nspins)
    println("inside neuromelanin nspins:   ", result.inside_neuromelanin_nspins)

    return result
end

# ============================================================
# Main script
# ============================================================

function main()
    println("============================================================")
    println("02 — Run cell + neuromelanin gradient-echo simulation")
    println("============================================================")

    ensure_output_dirs()

    println("\nConfiguration")
    println("-------------")
    println("cell SWC:              ", CELL_SWC_NAME)
    println("neuromelanin SWC:      ", NM_SWC_NAME)
    println("output prefix:         ", OUTPUT_PREFIX)
    println("nspins:                ", NSPINS)
    println("bbox margin:           ", BBOX_MARGIN)
    println("random seed:           ", RANDOM_SEED)
    println("TE:                    ", TE, " ms")
    println("TR:                    ", TR, " ms")
    println("diffusivity:           ", DIFFUSIVITY)
    println("skip_TR:               ", SKIP_TR)
    println("global R1:             ", GLOBAL_R1, " ms^-1")
    println("global R2:             ", GLOBAL_R2, " ms^-1")
    println("cell R1_inside:        ", CELL_R1_INSIDE, " ms^-1")
    println("NM R1 values:          ", NM_R1_INSIDE_VALUES, " ms^-1")
    println("cell permeability:     ", CELL_PERMEABILITY)
    println("NM permeability:       ", NM_PERMEABILITY)

    println("\nSelected input files")
    println("--------------------")
    print_selected_substrates(
        CELL_SWC_NAME,
        NM_SWC_NAME,
    )

    output_csv = gradient_echo_signal_output_path(OUTPUT_PREFIX)

    rows = Any[]
    push!(rows, gradient_echo_signal_header())

    for nm_R1_inside in NM_R1_INSIDE_VALUES
        result = run_one_nm_R1_condition(nm_R1_inside)

        append_gradient_echo_signal_row!(
            rows,
            CELL_SWC_NAME,
            NM_SWC_NAME,
            TE,
            TR,
            DIFFUSIVITY,
            GLOBAL_R1,
            GLOBAL_R2,
            CELL_R1_INSIDE,
            nm_R1_inside,
            NSPINS,
            SKIP_TR,
            result,
        )

        # Save progressively after each condition.
        write_rows_csv(
            output_csv,
            rows,
        )
    end

    println("\n============================================================")
    println("Simulation completed successfully.")
    println("============================================================")

    println("\nOutput CSV:")
    println("  ", output_csv)

    println("\nNext step:")
    println("  Run the plotting script:")
    println("  julia +1.11 --project=. docs/scripts/06_two_substrate_cell_neuromelanin/03_plot_cell_neuromelanin_signals.jl")
end

main()
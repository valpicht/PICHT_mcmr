# ============================================================
# 01_run_single_substrate_DWI.jl
#
# Goal:
# Run a reusable single-substrate DWI simulation.
#
# The substrate can be provided as:
#   - SWC file: id type x y z radius parent
#   - CSV file: x,y,z,radius
#
# The substrate is converted into an overlapping-sphere MCMR geometry.
#
# This script:
#   1. loads the substrate;
#   2. creates an overlapping-sphere geometry;
#   3. generates random spins in a bounding box around the substrate;
#   4. optionally keeps only inside spins;
#   5. runs DWI simulations for several TE values and directions;
#   6. computes ln(S/S0) and ADC;
#   7. saves detailed signal and ADC summary CSV files.
#
# Run from the repository root with:
#
#   julia +1.11 --project=. docs/scripts/05_general_substrate_simulation/01_run_single_substrate_dwi.jl
# ============================================================

using MCMRSimulator
using MRIBuilder
using DelimitedFiles
using Random
using Statistics

include(joinpath(@__DIR__, "utils_general_simulation.jl"))

# ============================================================
# User parameters
# ============================================================

# Substrate file name, without extension.
#
# The file must be located in:
#
#   substrate_generation/example_outputs/
#
# Examples:
#   SUBSTRATE_NAME = "one_soma_list_new"
#   SUBSTRATE_NAME = "cell_soma_r19_one_branche_neurons_list_new"
#   SUBSTRATE_NAME = "one_soma_with_2_branches_list_new"
#
const SUBSTRATE_NAME = "cell_soma_r19_one_branche_neurons_list_new"

# Supported formats:
#   "swc"
#   "csv"
const SUBSTRATE_FORMAT = "swc"

# Simulation output prefix.
#
# The script will create:
#
#   docs/data/processed/<OUTPUT_PREFIX>_signal.csv
#   docs/data/processed/<OUTPUT_PREFIX>_ADC_summary.csv
#
const OUTPUT_PREFIX = "single_substrate_DWI"

# Diffusion-weighting parameters.
const BVALS = [0.0, 0.5, 1.0]

# Echo times to test.
const TE_VALUES = [40]

# Repetition time.
const TR = 150

# Diffusivity in um^2/ms.
const DIFFUSIVITY = 1.0

# Number of spins generated in the substrate bounding box.
const NSPINS = 1_000

# Number of independent repeats.
#
# Each repeat uses a different random seed.
const N_REPEATS = 1

# Bounding box margin around the substrate.
const BBOX_MARGIN = 2.0

# Base random seed.
const RANDOM_SEED = 1234

# Manual timestep.
#
# If set to nothing, MCMR chooses the timestep automatically.
# If set to a number, this value is passed to Simulation(..., timestep = ...).
#
# A manual timestep can make simulations more feasible for very small spheres,
# but it is an approximation and should be reported.
const MANUAL_TIMESTEP = 0.01  # ms
# const MANUAL_TIMESTEP = nothing

# Whether to simulate only spins inside the substrate.
#
# true:
#   Uses only inside spins. Useful for directional diffusion validation.
#
# false:
#   Uses all spins in the bounding box.
#
const USE_INSIDE_SPINS_ONLY = true

# Diffusion directions to test.
#
# For an isotropic or simple substrate, one direction may be enough.
# For directional validation, use x/y/z.
const DIFFUSION_DIRECTIONS = [
    ("x", [1.0, 0.0, 0.0]),
    ("y", [0.0, 1.0, 0.0]),
    ("z", [0.0, 0.0, 1.0]),
]

# ============================================================
# Main script
# ============================================================

function main()
    println("============================================================")
    println("01 — Run single-substrate DWI simulation")
    println("============================================================")

    println("\nConfiguration")
    println("-------------")
    println("substrate name:   ", SUBSTRATE_NAME)
    println("substrate format: ", SUBSTRATE_FORMAT)
    println("output prefix:    ", OUTPUT_PREFIX)
    println("b-values:         ", BVALS)
    println("TE values:        ", TE_VALUES)
    println("TR:               ", TR, " ms")
    println("diffusivity:      ", DIFFUSIVITY)
    println("nspins:           ", NSPINS)
    println("n repeats:        ", N_REPEATS)
    println("bbox margin:      ", BBOX_MARGIN)
    println("random seed:      ", RANDOM_SEED)
    println("manual timestep:  ", MANUAL_TIMESTEP)
    println("inside only:      ", USE_INSIDE_SPINS_ONLY)

    println("\nDiffusion directions:")
    for (direction_name, direction) in DIFFUSION_DIRECTIONS
        println("  ", direction_name, " -> ", direction)
    end

    print_selected_substrate(
        SUBSTRATE_NAME,
        SUBSTRATE_FORMAT,
    )

    signal_csv, summary_csv = make_output_paths(OUTPUT_PREFIX)

    signal_rows = Any[]
    summary_rows = Any[]

    push!(signal_rows, signal_header())
    push!(summary_rows, summary_header())

    println("\nLoading substrate and creating overlapping-sphere geometry...")
    geometry, positions, radii, substrate_path = load_overlapping_substrate(
        SUBSTRATE_NAME,
        SUBSTRATE_FORMAT,
    )

    println("\nSubstrate extent:")
    print_sphere_ranges(positions, radii)

    for repeat_index in 1:N_REPEATS
        println("\n============================================================")
        println("Repeat ", repeat_index, " / ", N_REPEATS)
        println("============================================================")

        repeat_seed = RANDOM_SEED + repeat_index - 1

        println("Repeat seed: ", repeat_seed)

        all_spins, bbox = create_spins(
            positions,
            radii;
            nspins = NSPINS,
            margin = BBOX_MARGIN,
            seed = repeat_seed,
        )

        n_inside, n_outside = print_inside_outside_counts(
            geometry,
            all_spins,
        )

        inside_spins, outside_spins = split_inside_outside_spins(
            geometry,
            all_spins,
        )

        if USE_INSIDE_SPINS_ONLY
            if length(inside_spins) == 0
                error("""
                No inside spins were generated.

                Try increasing NSPINS or BBOX_MARGIN, or check the substrate geometry.
                """)
            end

            spins_to_simulate = inside_spins

            println("Using inside spins only.")
            println("Number of simulated spins: ", length(spins_to_simulate))
        else
            spins_to_simulate = all_spins

            println("Using all generated spins.")
            println("Number of simulated spins: ", length(spins_to_simulate))
        end

        for TE in TE_VALUES
            println("\n------------------------------------------------------------")
            println("Processing TE = ", TE, " ms")
            println("------------------------------------------------------------")

            for (direction_name, direction) in DIFFUSION_DIRECTIONS
                signals, ln_signal, adc, intercept, fit_ln =
                    run_single_TE_direction(
                        geometry,
                        spins_to_simulate,
                        BVALS,
                        TE,
                        direction_name,
                        direction;
                        TR = TR,
                        diffusivity = DIFFUSIVITY,
                        timestep = MANUAL_TIMESTEP,
                        skip_TR = 0,
                    )

                append_signal_rows!(
                    signal_rows,
                    SUBSTRATE_NAME,
                    SUBSTRATE_FORMAT,
                    TE,
                    direction_name,
                    direction,
                    BVALS,
                    repeat_index,
                    length(spins_to_simulate),
                    n_inside,
                    n_outside,
                    signals,
                    ln_signal,
                    fit_ln,
                )

                append_summary_row!(
                    summary_rows,
                    SUBSTRATE_NAME,
                    SUBSTRATE_FORMAT,
                    TE,
                    direction_name,
                    direction,
                    repeat_index,
                    length(spins_to_simulate),
                    n_inside,
                    n_outside,
                    adc,
                    intercept,
                )

                # Save progressively, so results are not lost if a later
                # direction or TE takes too long.
                write_rows_csv(signal_csv, signal_rows)
                write_rows_csv(summary_csv, summary_rows)
            end
        end
    end

    println("\n============================================================")
    println("Simulation completed successfully.")
    println("============================================================")

    println("\nOutputs:")
    println("  Signal CSV:  ", signal_csv)
    println("  Summary CSV: ", summary_csv)

    println("\nNext step:")
    println("  Run the plotting script:")
    println("  julia +1.11 --project=. docs/scripts/05_general_substrate_simulation/02_plot_single_substrate_dwi.jl")
end

main()
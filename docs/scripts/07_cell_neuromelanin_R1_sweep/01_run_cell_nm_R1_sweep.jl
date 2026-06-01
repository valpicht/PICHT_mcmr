# ============================================================
# 01_run_cell_nm_R1_sweep.jl
#
# Goal:
# Run a cell + neuromelanin GradientEcho simulation while sweeping
# the total R1 assigned to the neuromelanin compartment.
#
# This script follows the supervisor's request:
#
#   R1_NM = 1, 10, 20, 30, ..., 100, 200, 300, ..., 1000 s^-1
#
# Then it computes the NM-MRI contrast ratio:
#
#   CR = (S_i - S_0) / S_0
#
# where:
#   S_i = signal for a given R1_NM value
#   S_0 = baseline signal when R1_NM is equal to the background R1
#
# The script saves:
#
#   docs/data/processed/cell_nm_R1_sweep_R1_sweep.csv
#   docs/data/processed/cell_nm_R1_sweep_local_R1_check.csv
#
# Run from the repository root with:
#
#   julia +1.11 --project=. docs/scripts/07_cell_neuromelanin_R1_sweep/01_run_cell_nm_R1_sweep.jl
# ============================================================

using MCMRSimulator
using MRIBuilder
using DelimitedFiles
using Random
using Statistics

include(joinpath(@__DIR__, "utils_cell_nm_R1_sweep.jl"))

# ============================================================
# User parameters
# ============================================================

# The files must be located in:
#
#   substrate_generation/example_outputs/
#
# Example expected files:
#
#   substrate_generation/example_outputs/cell_soma_r19_one_branche_neurons_list_new.swc
#   substrate_generation/example_outputs/neuromelanin_sphere_r15_neurons_list_new.swc
#
# Change these names if your current files have different names.
const CELL_SWC_NAME = "cell_soma_r19_one_branche_neurons_list_new"
const NM_SWC_NAME = "neuromelanin_sphere_r15_neurons_list_new"

# Output prefix.
#
# The script will create:
#
#   docs/data/processed/<OUTPUT_PREFIX>_R1_sweep.csv
#   docs/data/processed/<OUTPUT_PREFIX>_local_R1_check.csv
#
const OUTPUT_PREFIX = "cell_nm_R1_sweep"

# Global/background R1.
#
# Unit: s^-1
#
# Baseline signal S0 is defined as the signal where:
#
#   R1_NM_total_s = R1_BACKGROUND_S
#
const R1_BACKGROUND_S = 1.0

# Global R2.
#
# Unit: s^-1
#
# If your previous simulation used "global R2 = 1/3 ms^-1",
# then that would correspond to 333.333 s^-1.
#
# In this script the parameter is in s^-1 and is converted internally
# to ms^-1.
#
# Check carefully which convention you used before.
const R2_GLOBAL_S = 1.0 / 3.0

# Neuromelanin R1 sweep.
#
# Full version suggested by the supervisor:
#
#   1, 10, 20, ..., 100, 200, 300, ..., 1000 s^-1
#
const R1_NM_TOTAL_VALUES_S = default_R1_NM_values_s()

# If the full version takes too long, use the coarser version instead:
#
# const R1_NM_TOTAL_VALUES_S = coarse_R1_NM_values_s()

# GradientEcho sequence parameters.
#
# Unit: ms
const TE = 3
const TR = 25

# Diffusion parameter.
const DIFFUSIVITY = 1.0

# Number of spins.
#
# Use a smaller value for debugging, then increase.
const NSPINS = 50_000

# Random seed for reproducibility.
const RANDOM_SEED = 1234

# Bounding box margin around the cell + NM substrates.
const BBOX_MARGIN = 2.0

# Number of TRs skipped before readout.
#
# Your previous gradient-echo tests often used skip_TR = 20.
const SKIP_TR = 20

# Manual timestep.
#
# Set to `nothing` to let MCMR choose automatically.
#
# If the simulation is too slow, you can test a manual timestep, but then
# report it clearly.
const MANUAL_TIMESTEP = nothing
# const MANUAL_TIMESTEP = 0.01

# Whether to write intermediate CSV files after each R1_NM value.
#
# Recommended: true, so results are saved even if a long sweep is interrupted.
const SAVE_PROGRESSIVELY = true

# ============================================================
# Helper functions
# ============================================================

"""
    ordered_R1_values(values, baseline)

Return R1 values with the baseline first.

This ensures that the baseline signal S0 is computed before all contrast
ratios are calculated.
"""
function ordered_R1_values(values, baseline)
    ordered = Float64[]

    push!(ordered, baseline)

    for value in values
        if value != baseline
            push!(ordered, value)
        end
    end

    return ordered
end

"""
    print_configuration()

Print all main simulation parameters.
"""
function print_configuration()
    println("Configuration")
    println("-------------")
    println("Cell SWC name:       ", CELL_SWC_NAME)
    println("NM SWC name:         ", NM_SWC_NAME)
    println("Output prefix:       ", OUTPUT_PREFIX)
    println("R1 background:       ", R1_BACKGROUND_S, " s^-1")
    println("R2 global:           ", R2_GLOBAL_S, " s^-1")
    println("TE:                  ", TE, " ms")
    println("TR:                  ", TR, " ms")
    println("diffusivity:         ", DIFFUSIVITY)
    println("nspins:              ", NSPINS)
    println("random seed:         ", RANDOM_SEED)
    println("bbox margin:         ", BBOX_MARGIN)
    println("skip_TR:             ", SKIP_TR)
    println("manual timestep:     ", MANUAL_TIMESTEP)
    println("save progressively:  ", SAVE_PROGRESSIVELY)
    println("")
    println("R1_NM total values:")
    println(R1_NM_TOTAL_VALUES_S)
end

"""
    build_check_positions(cell_positions, nm_positions)

Create a small list of points for the local R1 check.

The first point is the first cell sphere center.
The second point is the first NM sphere center.

This is a lightweight diagnostic to verify whether these points are inside
the expected geometries.
"""
function build_check_positions(cell_positions, nm_positions)
    check_positions = Vector{Vector{Float64}}()

    push!(check_positions, cell_positions[1])
    push!(check_positions, nm_positions[1])

    return check_positions
end

# ============================================================
# Main script
# ============================================================

function main()
    println("============================================================")
    println("01 — Run cell + neuromelanin R1 sweep")
    println("============================================================")

    print_configuration()

    print_selected_files(
        CELL_SWC_NAME,
        NM_SWC_NAME,
    )

    sweep_csv, local_R1_check_csv = make_output_paths(OUTPUT_PREFIX)

    sweep_rows = Any[]
    local_R1_check_rows = Any[]

    push!(sweep_rows, R1_sweep_header())
    push!(local_R1_check_rows, local_R1_check_header())

    println("\nLoading substrates")
    println("------------------")

    cell_positions, cell_radii, cell_swc_path = load_swc_spheres(
        CELL_SWC_NAME,
    )

    nm_positions, nm_radii, nm_swc_path = load_swc_spheres(
        NM_SWC_NAME,
    )

    println("\nSubstrate spatial extents")
    println("-------------------------")

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

    println("\nCreating common spin distribution")
    println("---------------------------------")

    spins, bbox = create_spins_for_cell_nm(
        cell_positions,
        cell_radii,
        nm_positions,
        nm_radii;
        nspins = NSPINS,
        margin = BBOX_MARGIN,
        seed = RANDOM_SEED,
    )

    # Create a baseline geometry only for initial compartment counting.
    baseline_geometry, baseline_R1_additional_s, baseline_R1_additional_ms =
        create_cell_nm_geometry(
            cell_positions,
            cell_radii,
            nm_positions,
            nm_radii;
            R1_background_s = R1_BACKGROUND_S,
            R1_NM_total_s = R1_BACKGROUND_S,
        )

    baseline_cell_geometry = baseline_geometry[1]
    baseline_nm_geometry = baseline_geometry[2]

    println("\nInitial compartment counts")
    println("--------------------------")

    print_compartment_counts(
        baseline_cell_geometry,
        baseline_nm_geometry,
        spins,
    )

    check_positions = build_check_positions(
        cell_positions,
        nm_positions,
    )

    println("\nR1 sweep")
    println("--------")

    R1_values_ordered = ordered_R1_values(
        R1_NM_TOTAL_VALUES_S,
        R1_BACKGROUND_S,
    )

    baseline_signals = nothing

    for (index, R1_NM_total_s) in enumerate(R1_values_ordered)
        println("\n============================================================")
        println("R1_NM value ", index, " / ", length(R1_values_ordered))
        println("R1_NM_total_s = ", R1_NM_total_s, " s^-1")
        println("============================================================")

        signals, R1_NM_additional_s, R1_NM_additional_ms, geometry =
            @time run_one_R1_NM_simulation(
                spins,
                cell_positions,
                cell_radii,
                nm_positions,
                nm_radii;
                R1_background_s = R1_BACKGROUND_S,
                R1_NM_total_s = R1_NM_total_s,
                R2_global_s = R2_GLOBAL_S,
                TE = TE,
                TR = TR,
                diffusivity = DIFFUSIVITY,
                timestep = MANUAL_TIMESTEP,
                skip_TR = SKIP_TR,
            )

        cell_geometry = geometry[1]
        nm_geometry = geometry[2]

        println("\nR1 values")
        println("---------")
        println("R1 background:       ", R1_BACKGROUND_S, " s^-1")
        println("R1 NM total:         ", R1_NM_total_s, " s^-1")
        println("R1 NM additional:    ", R1_NM_additional_s, " s^-1")
        println("R1 NM additional:    ", R1_NM_additional_ms, " ms^-1")

        println("\nSignals")
        println("-------")
        println("total:        ", signals.total)
        println("inside cell:  ", signals.inside_cell)
        println("outside cell: ", signals.outside_cell)
        println("inside NM:    ", signals.inside_nm)
        println("outside NM:   ", signals.outside_nm)

        if R1_NM_total_s == R1_BACKGROUND_S
            baseline_signals = signals

            println("\nBaseline S0 set from R1_NM_total_s = ", R1_NM_total_s, " s^-1")
        end

        if isnothing(baseline_signals)
            error("""
            Baseline signals have not been computed yet.

            The R1 sweep must include R1_BACKGROUND_S as the baseline.
            Current R1_BACKGROUND_S:
                $R1_BACKGROUND_S
            """)
        end

        CRs = compute_CRs(
            signals,
            baseline_signals,
        )

        println("\nContrast ratios relative to baseline")
        println("------------------------------------")
        println("CR total:        ", CRs.CR_total)
        println("CR inside cell:  ", CRs.CR_inside_cell)
        println("CR outside cell: ", CRs.CR_outside_cell)
        println("CR inside NM:    ", CRs.CR_inside_nm)
        println("CR outside NM:   ", CRs.CR_outside_nm)

        append_R1_sweep_row!(
            sweep_rows,
            CELL_SWC_NAME,
            NM_SWC_NAME,
            R1_BACKGROUND_S,
            R1_NM_total_s,
            R1_NM_additional_s,
            R1_NM_additional_ms,
            R2_GLOBAL_S,
            TE,
            TR,
            DIFFUSIVITY,
            MANUAL_TIMESTEP,
            SKIP_TR,
            signals,
            CRs,
        )

        append_local_R1_check_rows!(
            local_R1_check_rows,
            CELL_SWC_NAME,
            NM_SWC_NAME,
            R1_BACKGROUND_S,
            R1_NM_total_s,
            R1_NM_additional_s,
            cell_geometry,
            nm_geometry,
            check_positions,
        )

        if SAVE_PROGRESSIVELY
            write_rows_csv(
                sweep_csv,
                sweep_rows,
            )

            write_rows_csv(
                local_R1_check_csv,
                local_R1_check_rows,
            )

            println("Progress saved.")
        end
    end

    if !SAVE_PROGRESSIVELY
        write_rows_csv(
            sweep_csv,
            sweep_rows,
        )

        write_rows_csv(
            local_R1_check_csv,
            local_R1_check_rows,
        )
    end

    println("\n============================================================")
    println("R1 sweep completed successfully.")
    println("============================================================")

    println("\nOutput files:")
    println("  Sweep CSV:")
    println("    ", sweep_csv)
    println("  Local R1 check CSV:")
    println("    ", local_R1_check_csv)

    println("\nNext step:")
    println("  Plot the contrast ratio:")
    println("  julia +1.11 --project=. docs/scripts/07_cell_neuromelanin_R1_sweep/02_plot_cell_nm_R1_sweep_CR.jl")
end

main()
# ============================================================
# 02_run_dwi_overlapping_spheres.jl
#
# Goal:
# Run a first diffusion-weighted simulation on an overlapping-sphere
# substrate generated from an SWC file.
#
# This script:
#   1. loads an SWC substrate as overlapping spheres;
#   2. creates random spins around the substrate;
#   3. runs DWI simulations for several b-values;
#   4. separates inside and outside signals;
#   5. saves raw and normalized signals to CSV.
#
# Run from the repository root with:
#
#     julia +1.11 --project=. docs/scripts/04_overlapping_spheres/02_run_dwi_overlapping_spheres.jl
# ============================================================

using MCMRSimulator
using MRIBuilder
using CairoMakie
using DelimitedFiles

include(joinpath(@__DIR__, "utils_overlapping.jl"))

# ============================================================
# Parameters
# ============================================================

const SWC_NAME = "one_soma_list_new"

const NSPINS = 10_000
const BBOX_MARGIN = 2.0

const BVALS = collect(0.0:0.1:1.0)
const TE = 80
const TR = 300
const DIFFUSIVITY = 1.0
const SKIP_TR = 2

# ============================================================
# Main script
# ============================================================

function main()
    println("============================================================")
    println("02 — Run DWI on overlapping spheres")
    println("============================================================")

    output_dir, figure_dir = ensure_output_dirs()

    output_csv = joinpath(
        output_dir,
        "overlapping_spheres_signal_vs_b_data.csv",
    )

    inside_fig_path = joinpath(
        figure_dir,
        "overlapping_spheres_inside_signal_vs_b.png",
    )

    outside_fig_path = joinpath(
        figure_dir,
        "overlapping_spheres_outside_signal_vs_b.png",
    )

    println("\nSelected SWC substrate:")
    println(SWC_NAME)

    geometry, positions, radii = load_overlapping_spheres(SWC_NAME)

    println("\nSphere coordinate ranges:")
    print_sphere_ranges(positions, radii)

    println("\nCreating spins...")
    spins = create_spins_for_spheres(
        positions,
        radii;
        nspins = NSPINS,
        margin = BBOX_MARGIN,
    )

    println("\nInside/outside spin count before simulation:")
    print_inside_outside_counts(geometry, spins)

    println("\nSimulation parameters")
    println("---------------------")
    println("b-values: ", BVALS)
    println("TE = ", TE, " ms")
    println("TR = ", TR, " ms")
    println("diffusivity = ", DIFFUSIVITY)
    println("skip_TR = ", SKIP_TR)
    println("nspins = ", NSPINS)

    println("\nRunning b-value sweep...")
    inside_signals, outside_signals = simulate_signal_vs_b(
        spins,
        geometry,
        BVALS;
        TE = TE,
        TR = TR,
        diffusivity = DIFFUSIVITY,
        skip_TR = SKIP_TR,
    )

    inside_norm = normalize_by_first_value(inside_signals)
    outside_norm = normalize_by_first_value(outside_signals)

    println("\nInside normalized signal:")
    for (b, s) in zip(BVALS, inside_norm)
        println("b = ", b, " -> ", s)
    end

    println("\nOutside normalized signal:")
    for (b, s) in zip(BVALS, outside_norm)
        println("b = ", b, " -> ", s)
    end

    println("\nSaving CSV output...")
    write_signal_vs_b_csv(
        output_csv,
        BVALS,
        inside_signals,
        outside_signals,
    )

    println("\nPlotting inside signal...")
    fig_inside = Figure()
    ax_inside = Axis(
        fig_inside[1, 1],
        title = "Overlapping spheres: inside signal",
        xlabel = "b",
        ylabel = "S(b)/S(0)",
    )

    lines!(ax_inside, BVALS, inside_norm)
    scatter!(ax_inside, BVALS, inside_norm)

    save(inside_fig_path, fig_inside)
    display(fig_inside)

    println("Inside signal figure saved:")
    println(inside_fig_path)

    println("\nPlotting outside signal...")
    fig_outside = Figure()
    ax_outside = Axis(
        fig_outside[1, 1],
        title = "Overlapping spheres: outside signal",
        xlabel = "b",
        ylabel = "S(b)/S(0)",
    )

    lines!(ax_outside, BVALS, outside_norm)
    scatter!(ax_outside, BVALS, outside_norm)

    save(outside_fig_path, fig_outside)
    display(fig_outside)

    println("Outside signal figure saved:")
    println(outside_fig_path)

    println("\nOverlapping-sphere DWI simulation completed successfully.")
end

main()
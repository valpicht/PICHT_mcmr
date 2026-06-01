# ============================================================
# 03_compare_non_overlapping_vs_overlapping.jl
#
# Goal:
# Compare the behavior of the same SWC-derived sphere substrate with
# and without the overlapping-sphere option.
#
# This script checks the effect of:
#
#     overlapping = false
#     overlapping = true
#
# for the same sphere positions, radii, spin distribution, and DWI
# sequence parameters.
#
# This is not yet a mesh-vs-overlapping comparison. That comparison is
# done later in:
#
#     docs/scripts/05_mesh_vs_overlapping_comparison/
#
# Run from the repository root with:
#
#     julia +1.11 --project=. docs/scripts/04_overlapping_spheres/03_compare_non_overlapping_vs_overlapping.jl
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
# Helper functions
# ============================================================

"""
    create_spheres_geometry(positions, radii; overlapping)

Create an MCMR Spheres geometry from positions and radii.

The `overlapping` argument controls whether overlapping spheres are
interpreted as a valid continuous substrate.
"""
function create_spheres_geometry(positions, radii; overlapping::Bool)
    return Spheres(
        radius = radii,
        position = positions,
        overlapping = overlapping,
    )
end

"""
    run_b_sweep_for_geometry(spins, geometry, label)

Run the same DWI b-value sweep for one geometry and return inside/outside
signals and normalized signals.
"""
function run_b_sweep_for_geometry(spins, geometry, label)
    println("\nRunning b-value sweep for: ", label)

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

    return inside_signals, outside_signals, inside_norm, outside_norm
end

# ============================================================
# Main script
# ============================================================

function main()
    println("============================================================")
    println("03 — Compare non-overlapping vs overlapping spheres")
    println("============================================================")

    output_dir, figure_dir = ensure_output_dirs()

    output_csv = joinpath(
        output_dir,
        "non_overlapping_vs_overlapping_signal_vs_b_data.csv",
    )

    inside_fig_path = joinpath(
        figure_dir,
        "non_overlapping_vs_overlapping_inside_signal_vs_b.png",
    )

    outside_fig_path = joinpath(
        figure_dir,
        "non_overlapping_vs_overlapping_outside_signal_vs_b.png",
    )

    swc_path = swc_path_from_name(SWC_NAME)

    println("\nSelected SWC substrate:")
    println(SWC_NAME)

    positions, radii = read_swc_spheres(swc_path)

    println("\nSphere coordinate ranges:")
    print_sphere_ranges(positions, radii)

    geometry_non_overlapping = create_spheres_geometry(
        positions,
        radii;
        overlapping = false,
    )

    geometry_overlapping = create_spheres_geometry(
        positions,
        radii;
        overlapping = true,
    )

    println("\nGeometry types:")
    println("Non-overlapping: ", typeof(geometry_non_overlapping))
    println("Overlapping:     ", typeof(geometry_overlapping))

    println("\nCreating common spin distribution...")
    spins = create_spins_for_spheres(
        positions,
        radii;
        nspins = NSPINS,
        margin = BBOX_MARGIN,
    )

    println("\nInside/outside count for non-overlapping geometry:")
    print_inside_outside_counts(geometry_non_overlapping, spins)

    println("\nInside/outside count for overlapping geometry:")
    print_inside_outside_counts(geometry_overlapping, spins)

    non_inside, non_outside, non_inside_norm, non_outside_norm =
        run_b_sweep_for_geometry(
            spins,
            geometry_non_overlapping,
            "overlapping = false",
        )

    ov_inside, ov_outside, ov_inside_norm, ov_outside_norm =
        run_b_sweep_for_geometry(
            spins,
            geometry_overlapping,
            "overlapping = true",
        )

    println("\nSaving CSV output...")

    data = hcat(
        BVALS,
        non_inside,
        non_inside_norm,
        non_outside,
        non_outside_norm,
        ov_inside,
        ov_inside_norm,
        ov_outside,
        ov_outside_norm,
    )

    open(output_csv, "w") do io
        writedlm(
            io,
            [
                "b" "nonoverlap_inside_signal" "nonoverlap_inside_norm" "nonoverlap_outside_signal" "nonoverlap_outside_norm" "overlap_inside_signal" "overlap_inside_norm" "overlap_outside_signal" "overlap_outside_norm"
            ],
            ',',
        )
        writedlm(io, data, ',')
    end

    println("CSV saved:")
    println(output_csv)

    println("\nPlotting inside signal comparison...")

    fig_inside = Figure()
    ax_inside = Axis(
        fig_inside[1, 1],
        title = "Inside signal: non-overlapping vs overlapping spheres",
        xlabel = "b",
        ylabel = "S(b)/S(0)",
    )

    lines!(
        ax_inside,
        BVALS,
        non_inside_norm,
        label = "overlapping = false",
    )
    scatter!(
        ax_inside,
        BVALS,
        non_inside_norm,
    )

    lines!(
        ax_inside,
        BVALS,
        ov_inside_norm,
        label = "overlapping = true",
    )
    scatter!(
        ax_inside,
        BVALS,
        ov_inside_norm,
    )

    axislegend(ax_inside, position = :rb)

    save(inside_fig_path, fig_inside)
    display(fig_inside)

    println("Inside comparison figure saved:")
    println(inside_fig_path)

    println("\nPlotting outside signal comparison...")

    fig_outside = Figure()
    ax_outside = Axis(
        fig_outside[1, 1],
        title = "Outside signal: non-overlapping vs overlapping spheres",
        xlabel = "b",
        ylabel = "S(b)/S(0)",
    )

    lines!(
        ax_outside,
        BVALS,
        non_outside_norm,
        label = "overlapping = false",
    )
    scatter!(
        ax_outside,
        BVALS,
        non_outside_norm,
    )

    lines!(
        ax_outside,
        BVALS,
        ov_outside_norm,
        label = "overlapping = true",
    )
    scatter!(
        ax_outside,
        BVALS,
        ov_outside_norm,
    )

    axislegend(ax_outside, position = :rb)

    save(outside_fig_path, fig_outside)
    display(fig_outside)

    println("Outside comparison figure saved:")
    println(outside_fig_path)

    println("\nSummary")
    println("-------")
    println("This script compares the same SWC-derived sphere substrate with")
    println("overlapping = false and overlapping = true.")
    println("For a substrate containing only one sphere, both modes are expected")
    println("to give similar or identical results. Differences become meaningful")
    println("for substrates containing multiple overlapping spheres.")

    println("\nNon-overlapping vs overlapping comparison completed successfully.")
end

main()
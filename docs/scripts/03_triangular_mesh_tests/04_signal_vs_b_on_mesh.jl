# ============================================================
# 04_signal_vs_b_on_mesh.jl
#
# Goal:
# Simulate the inside and outside signal of a triangular mesh
# substrate for several b-values.
#
# This script is a cleaned version of the original
# plot_signal_vs_b.jl script.
#
# It:
#   1. loads a triangular mesh;
#   2. creates random spins in a bounding box around the mesh;
#   3. separates spins into inside and outside compartments;
#   4. simulates the signal for several b-values;
#   5. normalizes the signal by S(0);
#   6. saves a CSV file and two figures.
#
# Run from the repository root with:
#
#     julia +1.11 --project=. docs/scripts/03_triangular_mesh_tests/04_signal_vs_b_on_mesh.jl
#
# Or run cell by cell in VS Code.
# ============================================================


# %% Cell 1 — Load packages and utilities

using MCMRSimulator
using MRIBuilder
using CairoMakie
using DelimitedFiles

include(joinpath(@__DIR__, "utils_mesh.jl"))


# %% Cell 2 — Choose mesh and output names

const MESH_NAME = "one_soma_list_new_mesh"

const OUTPUT_DIR = joinpath("docs", "data", "processed")
const FIGURE_DIR = joinpath("docs", "figures")

mkpath(OUTPUT_DIR)
mkpath(FIGURE_DIR)

const OUTPUT_CSV = joinpath(OUTPUT_DIR, "mesh_signal_vs_b_data.csv")
const INSIDE_FIG = joinpath(FIGURE_DIR, "mesh_inside_signal_vs_b.png")
const OUTSIDE_FIG = joinpath(FIGURE_DIR, "mesh_outside_signal_vs_b.png")

println("Selected mesh: ", MESH_NAME)


# %% Cell 3 — Define simulation helper

"""
    simulate_signal_vs_b(spins, geometry, bvals; TE = 40, TR = 300, diffusivity = 1.0)

Run one DWI simulation per b-value and return one scalar transverse
signal per b-value.

This function is intentionally close to the original project script:
for each b-value, it creates a DWI sequence, creates a Simulation, and
runs a readout after `skip_TR = 2`.
"""
function simulate_signal_vs_b(
    spins,
    geometry,
    bvals;
    TE = 40,
    TR = 300,
    diffusivity = 1.0,
)
    signals = Float64[]

    for b in bvals
        println("Running simulation for b = ", b)

        sequence = DWI(
            bval = b,
            TE = TE,
            TR = TR,
            scanner = MRIBuilder.Siemens_Prisma,
        )

        sim = Simulation(
            sequence,
            diffusivity = diffusivity,
            geometry = geometry,
        )

        sig = readout(
            spins,
            sim,
            skip_TR = 2,
        )

        # Original scripts used sig.orient.transverse.
        # abs(...) makes the scalar robust if the transverse signal is complex.
        push!(signals, abs(sig.orient.transverse))
    end

    return signals
end


# %% Cell 4 — Load mesh

geometry = load_example_mesh(MESH_NAME)

println("Mesh loaded.")


# %% Cell 5 — Create bounding box and spins

bbox = make_bbox_from_geometry(
    geometry;
    margin = 2.0,
)

println("Bounding box:")
println(bbox)

all_spins = Snapshot(10_000, bbox)

inside_mask = isinside(geometry, all_spins) .> 0
inside_spins = all_spins[inside_mask]
outside_spins = all_spins[.!inside_mask]

println("Number of inside spins: ", length(inside_spins))
println("Number of outside spins: ", length(outside_spins))


# %% Cell 6 — Define b-values

bvals = collect(0.0:0.1:1.0)

println("b-values:")
println(bvals)


# %% Cell 7 — Simulate inside and outside signal

inside_signals = simulate_signal_vs_b(
    inside_spins,
    geometry,
    bvals;
    TE = 40,
    TR = 300,
    diffusivity = 1.0,
)

outside_signals = simulate_signal_vs_b(
    outside_spins,
    geometry,
    bvals;
    TE = 40,
    TR = 300,
    diffusivity = 1.0,
)


# %% Cell 8 — Normalize signals by S(0)

inside_norm = inside_signals ./ inside_signals[1]
outside_norm = outside_signals ./ outside_signals[1]

println("\nInside normalized signal:")
for (b, s) in zip(bvals, inside_norm)
    println("b = ", b, " -> ", s)
end

println("\nOutside normalized signal:")
for (b, s) in zip(bvals, outside_norm)
    println("b = ", b, " -> ", s)
end


# %% Cell 9 — Save CSV output

data = hcat(
    bvals,
    inside_signals,
    inside_norm,
    outside_signals,
    outside_norm,
)

open(OUTPUT_CSV, "w") do io
    writedlm(
        io,
        ["b" "inside_signal" "inside_norm" "outside_signal" "outside_norm"],
        ',',
    )
    writedlm(io, data, ',')
end

println("\nData saved:")
println(" - ", OUTPUT_CSV)


# %% Cell 10 — Plot inside signal

fig_inside = Figure()
ax_inside = Axis(
    fig_inside[1, 1],
    title = "Inside mesh: S(b)/S(0) vs b",
    xlabel = "b",
    ylabel = "S(b)/S(0)",
)

lines!(ax_inside, bvals, inside_norm)
scatter!(ax_inside, bvals, inside_norm)

save(INSIDE_FIG, fig_inside)
display(fig_inside)

println("Inside signal figure saved:")
println(" - ", INSIDE_FIG)


# %% Cell 11 — Plot outside signal

fig_outside = Figure()
ax_outside = Axis(
    fig_outside[1, 1],
    title = "Outside mesh: S(b)/S(0) vs b",
    xlabel = "b",
    ylabel = "S(b)/S(0)",
)

lines!(ax_outside, bvals, outside_norm)
scatter!(ax_outside, bvals, outside_norm)

save(OUTSIDE_FIG, fig_outside)
display(fig_outside)

println("Outside signal figure saved:")
println(" - ", OUTSIDE_FIG)


# %% Cell 12 — End

println("\nSignal-vs-b mesh simulation completed successfully.")
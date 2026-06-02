# ============================================================
# 02_mesh_bounding_box_and_subsets.jl
#
# Goal:
# Inspect mesh vertices, create a bounding box, generate spins,
# and count how many spins are inside the mesh.
#
# This script is based on the original bounding-box debugging script.
#
# Run from the repository root with:
#
#     julia +1.11 --project=. docs/scripts/03_triangular_mesh_tests/02_mesh_bounding_box_and_subsets.jl
# ============================================================


# %% Cell 1 — Load packages

using MCMRSimulator
using MRIBuilder

include(joinpath(@__DIR__, "utils_mesh.jl"))


# %% Cell 2 — Choose and load mesh

const MESH_NAME = "one_soma_list_new_mesh"

geometry = load_example_mesh(MESH_NAME)


# %% Cell 3 — Inspect vertices

print_mesh_ranges(geometry)


# %% Cell 4 — Create bounding box

margin = 2.0

bbox = make_bbox_from_geometry(
    geometry;
    margin = margin,
)

println("Bounding box:")
println(bbox)


# %% Cell 5 — Generate spins

nspins = 100_000

all_spins = Snapshot(nspins, bbox)

println("Number of generated spins: ", length(all_spins))


# %% Cell 6 — Select spins inside the mesh

inside_spins = all_spins[isinside(geometry, all_spins) .> 0]

println("Number of spins inside geometry: ", length(inside_spins))


# %% Cell 7 — Create a simple DWI simulation

sequence = create_default_dwi_sequence(
    bval = 2.0,
    TE = 80,
    TR = 300,
)

sim = Simulation(
    sequence,
    diffusivity = 2.0,
    geometry = geometry,
)

println("Simulation created.")


# %% Cell 8 — Readout using all spins and inside/outside subsets

sig_comp = readout(
    all_spins,
    sim,
    skip_TR = 2,
    subset = [
        Subset(inside = true),
        Subset(inside = false),
    ],
)

println("Inside: ", sig_comp[1].nspins)
println("Outside: ", sig_comp[2].nspins)

println(sig_comp)


# %% Cell 9 — End

println("Mesh bounding box and subset test completed successfully.")
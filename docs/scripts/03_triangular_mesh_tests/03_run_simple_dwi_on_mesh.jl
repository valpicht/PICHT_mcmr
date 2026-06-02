# ============================================================
# 03_run_simple_dwi_on_mesh.jl
#
# Goal:
# Run a simple DWI simulation on a triangular mesh and separate
# inside/outside signal contributions.
#
# This is a cleaned version of the original no-branches mesh test.
#
# Run from the repository root with:
#
#     julia +1.11 --project=. docs/scripts/03_triangular_mesh_tests/03_run_simple_dwi_on_mesh.jl
# ============================================================


# %% Cell 1 — Load packages

using MCMRSimulator
using MRIBuilder

include(joinpath(@__DIR__, "utils_mesh.jl"))


# %% Cell 2 — Choose and load mesh

const MESH_NAME = "one_soma_list_new_mesh"

geometry = load_example_mesh(MESH_NAME)


# %% Cell 3 — Create DWI sequence and simulation

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


# %% Cell 4 — Super quick readout

sig_comp_small = readout(
    50,
    sim,
    skip_TR = 2,
    subset = [
        Subset(inside = true),
        Subset(inside = false),
    ],
)

println("SMALL Inside: ", sig_comp_small[1].nspins)
println("SMALL Outside: ", sig_comp_small[2].nspins)


# %% Cell 5 — Larger readout

sig_comp = readout(
    2000,
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


# %% Cell 6 — End

println("Simple DWI simulation on mesh completed successfully.")
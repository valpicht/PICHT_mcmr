# ============================================================
# 01_load_and_plot_mesh.jl
#
# Goal:
# Load a triangular mesh and plot it.
#
# This script is based on the first mesh-loading tests performed
# during the project, but uses a relative mesh path instead of an
# absolute path.
#
# Run from the repository root with:
#
#     julia +1.11 --project=. docs/scripts/03_triangular_mesh_tests/01_load_and_plot_mesh.jl
# ============================================================


# %% Cell 1 — Load packages

using MCMRSimulator
using CairoMakie

include(joinpath(@__DIR__, "utils_mesh.jl"))


# %% Cell 2 — Choose mesh

const MESH_NAME = "one_soma_list_new_mesh"

println("Selected mesh: ", MESH_NAME)


# %% Cell 3 — Load mesh

geometry = load_example_mesh(MESH_NAME)


# %% Cell 4 — Inspect mesh

print_mesh_ranges(geometry)


# %% Cell 5 — Plot mesh

# This keeps the same plotting style as the original working script.
f = plot(geometry)
display(f)


# %% Cell 6 — End

println("Mesh loading and plotting completed successfully.")
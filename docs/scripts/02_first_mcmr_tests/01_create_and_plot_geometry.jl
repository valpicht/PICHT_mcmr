# ============================================================
# 01_create_and_plot_geometry.jl
#
# Goal:
# First contact with MCMRSimulator.
#
# This script creates a simple cylinder geometry and plots it.
# It follows the first part of the official MCMR Julia tutorial.
#
# Run from the repository root with:
#
#     julia --project=. docs/scripts/02_first_mcmr_tests/01_create_and_plot_geometry.jl
#
# Or run cell by cell in VS Code.
# ============================================================


# %% Cell 1 — Load packages

using MCMRSimulator
using CairoMakie

# CairoMakie is used for static plots.
# GLMakie or WGLMakie can also be used for interactive plots.

update_theme!(Theme(Axis = (xgridvisible = false, ygridvisible = false)))


# %% Cell 2 — Create a simple geometry

geometry = Cylinders(
    radius = 1.0,
    repeats = [2.5, 2.5],
)

println("Geometry created.")
println("Geometry type: ", typeof(geometry))


# %% Cell 3 — Plot the geometry

plot_plane = PlotPlane(size = 5)

fig = plot(plot_plane, geometry)

xlims!(fig.axis, -2.5, 2.5)
ylims!(fig.axis, -2.5, 2.5)

display(fig)

println("Geometry plotted successfully.")


# %% Cell 4 — Why this geometry?

println("""
This geometry represents regularly packed cylinders.

In the MCMR tutorial, this is used as a simple model of packed axons.
The cylinder radius is 1 micrometer and the cylinders repeat every
2.5 micrometers in the x- and y-directions.
""")
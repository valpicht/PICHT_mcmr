# ============================================================
# 02_create_sequence_and_simulation.jl
#
# Goal:
# Create a diffusion-weighted sequence with MRIBuilder and combine
# it with a simple MCMR geometry into a Simulation object.
#
# Run from the repository root with:
#
#     julia --project=. docs/scripts/02_first_mcmr_tests/02_create_sequence_and_simulation.jl
#
# Or run cell by cell in VS Code.
# ============================================================


# %% Cell 1 — Load packages

using MCMRSimulator
using MRIBuilder
using CairoMakie

update_theme!(Theme(Axis = (xgridvisible = false, ygridvisible = false)))


# %% Cell 2 — Define geometry

geometry = Cylinders(
    radius = 1.0,
    repeats = [2.5, 2.5],
)

println("Geometry created.")
println("Geometry type: ", typeof(geometry))


# %% Cell 3 — Define a diffusion-weighted MRI sequence

sequence = DWI(
    bval = 2.0,
    TE = 80,
    TR = 300,
    scanner = Siemens_Prisma,
)

println("DWI sequence created.")
println("Sequence type: ", typeof(sequence))


# %% Cell 4 — Plot the sequence

fig_sequence = plot_sequence(sequence)
display(fig_sequence)

println("Sequence plotted successfully.")


# %% Cell 5 — Create a Simulation object

simulation = Simulation(
    sequence,
    R2 = 0.012,
    R1 = 3e-3,
    diffusivity = 2.0,
    off_resonance = 0.1,
    geometry = geometry,
)

println("Simulation created.")
println("Simulation type: ", typeof(simulation))


# %% Cell 6 — Explain the parameters

println("""
Simulation parameters used here:

- bval = 2.0
- TE = 80 ms
- TR = 300 ms
- diffusivity = 2.0 um^2/ms
- R1 = 3e-3 ms^-1
- R2 = 0.012 ms^-1
- off_resonance = 0.1

This is still a tutorial example. Later scripts modify these parameters
to study specific substrates and local relaxation effects.
""")
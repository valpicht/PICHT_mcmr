# ============================================================
# 04_signal_evolution_and_snapshots.jl
#
# Goal:
# Visualize signal evolution over multiple repetition times and return
# snapshots of spin states.
#
# This script follows the later parts of the official MCMR tutorial.
#
# Run from the repository root with:
#
#     julia --project=. docs/scripts/02_first_mcmr_tests/04_signal_evolution_and_snapshots.jl
#
# Or run cell by cell in VS Code.
# ============================================================


# %% Cell 1 — Load packages

using MCMRSimulator
using MRIBuilder
using CairoMakie

update_theme!(Theme(Axis = (xgridvisible = false, ygridvisible = false)))


# %% Cell 2 — Define geometry, sequence, and simulation

geometry = Cylinders(
    radius = 1.0,
    repeats = [2.5, 2.5],
)

sequence = DWI(
    bval = 2.0,
    TE = 80,
    TR = 300,
    scanner = Siemens_Prisma,
)

simulation = Simulation(
    sequence,
    R2 = 0.012,
    R1 = 3e-3,
    diffusivity = 2.0,
    off_resonance = 0.1,
    geometry = geometry,
)

println("Simulation ready.")


# %% Cell 3 — Signal over multiple TRs

signals = readout(
    1000,
    simulation,
    nTR = 6,
)

println("Number of readouts: ", length(signals))
println("Signals over TRs:")
println(signals)


# %% Cell 4 — Plot longitudinal and transverse components over TRs

fig = Figure()
ax = Axis(
    fig[1, 1],
    xlabel = "TR index",
    ylabel = "Signal",
    title = "Signal evolution over repeated TRs",
)

lines!(ax, longitudinal.(signals), label = "Longitudinal")
lines!(ax, transverse.(signals), label = "Transverse")
axislegend(ax)

display(fig)


# %% Cell 5 — Signal sampled at custom time points

times = 0:0.1:100

average_signals = readout(
    3000,
    simulation,
    times,
)

transverse_signal = transverse.(average_signals) ./ 3000

fig_time = Figure()
ax_time = Axis(
    fig_time[1, 1],
    xlabel = "Time (ms)",
    ylabel = "Mean transverse signal",
    title = "Transverse signal evolution",
)

lines!(ax_time, times, transverse_signal)

display(fig_time)


# %% Cell 6 — Return trajectory snapshots for a few spins

snapshots = readout(
    [[0, 0, 0], [1, 1, 0]],
    simulation,
    0:0.01:3.0,
    return_snapshot = true,
)

println("Trajectory snapshots created.")
println("Snapshot object type: ", typeof(snapshots))


# %% Cell 7 — Plot spin trajectories in 3D

fig_snapshots = plot(snapshots)
display(fig_snapshots)


# %% Cell 8 — Return a full snapshot at readout time

snapshot = readout(
    3000,
    simulation,
    return_snapshot = true,
)

plot_plane = PlotPlane(size = 2.5)

fig_snapshot = plot(plot_plane, snapshot)
plot!(plot_plane, geometry)

display(fig_snapshot)


# %% Cell 9 — Alternative snapshot visualization: dyads

fig_dyad = plot(plot_plane, snapshot, kind = :dyad)
plot!(plot_plane, geometry)

display(fig_dyad)


# %% Cell 10 — Alternative snapshot visualization: image

fig_image = plot(plot_plane, snapshot, kind = :image)
plot!(plot_plane, geometry)

display(fig_image)


# %% Cell 11 — End of script

println("Signal evolution and snapshots completed successfully.")
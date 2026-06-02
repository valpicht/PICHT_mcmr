# ============================================================
# 03_simple_readouts_and_subsets.jl
#
# Goal:
# Run simple MCMR readouts and split the signal into inside/outside
# compartments.
#
# This script introduces:
#   - readout(nspins, simulation)
#   - skip_TR
#   - Subset(inside=true)
#   - Subset(inside=false)
#
# Run from the repository root with:
#
#     julia --project=. docs/scripts/02_first_mcmr_tests/03_simple_readouts_and_subsets.jl
#
# Or run cell by cell in VS Code.
# ============================================================


# %% Cell 1 — Load packages

using MCMRSimulator
using MRIBuilder


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


# %% Cell 3 — First simple readout

signal = readout(1000, simulation)

println("Signal without skipped TR:")
println(signal)

println("""
This is the simplest readout call. It simulates 1000 spins and returns
the signal at the readout time defined by the sequence.
""")


# %% Cell 4 — Readout after skipped TRs

signal_after_skip = readout(
    1000,
    simulation,
    skip_TR = 2,
)

println("Signal after skip_TR = 2:")
println(signal_after_skip)

println("""
skip_TR lets the system evolve for several repetition times before
measuring the signal. This is useful because the first repetition may
not yet represent a steady or equilibrated signal.
""")


# %% Cell 5 — Compartment-wise readout

signals_by_compartment = readout(
    10000,
    simulation,
    skip_TR = 2,
    subset = [
        Subset(inside = true),
        Subset(inside = false),
    ],
)

inside_signal = signals_by_compartment[1]
outside_signal = signals_by_compartment[2]

println("Inside-cylinder signal:")
println(inside_signal)

println("\nOutside-cylinder signal:")
println(outside_signal)

println("\nNumber of spins:")
println("Inside: ", inside_signal.nspins)
println("Outside: ", outside_signal.nspins)


# %% Cell 6 — Extract longitudinal and transverse components

println("Inside transverse signal: ", transverse(inside_signal))
println("Outside transverse signal: ", transverse(outside_signal))

println("Inside longitudinal signal: ", longitudinal(inside_signal))
println("Outside longitudinal signal: ", longitudinal(outside_signal))

println("""
The transverse component is usually the most relevant quantity for
MRI signal readouts. The longitudinal component describes the remaining
z-magnetization.
""")


# %% Cell 7 — End of script

println("Simple readouts and subsets completed successfully.")
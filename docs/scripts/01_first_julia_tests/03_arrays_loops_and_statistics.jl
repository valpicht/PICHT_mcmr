# ============================================================
# 03_arrays_loops_and_statistics.jl
#
# Goal:
# This script introduces arrays, loops, dictionaries, simple statistics,
# and basic file output.
#
# These concepts are used later to:
#   - loop over b-values and echo times
#   - store simulation outputs
#   - compute simple summaries
#   - save results to CSV files
#
# Run from the repository root with:
#
#     julia --project=. docs/scripts/01_first_julia_tests/03_arrays_loops_and_statistics.jl
#
# Or run cell by cell in VS Code.
# ============================================================


# %% Cell 1 — Import useful standard packages

using Statistics
using DelimitedFiles

# Statistics provides mean(), std(), etc.
# DelimitedFiles provides writedlm(), useful for simple CSV-like files.


# %% Cell 2 — Vectors

bvals = [0.0, 0.5, 1.0, 1.5, 2.0]
signals = [1.0, 0.82, 0.67, 0.55, 0.45]

println("b-values:")
println(bvals)

println("\nSignals:")
println(signals)


# %% Cell 3 — Looping over one vector

println("Looping over b-values:")

for b in bvals
    println("Current b-value = ", b)
end


# %% Cell 4 — Looping over two vectors together

println("Looping over b-values and signals:")

for (b, signal) in zip(bvals, signals)
    println("b = ", b, " | signal = ", signal)
end

# zip() is useful when two arrays represent matching values.


# %% Cell 5 — Preallocating an array

normalized_signals = zeros(length(signals))

for i in eachindex(signals)
    normalized_signals[i] = signals[i] / signals[1]
end

println("Normalized signals:")
println(normalized_signals)

# Preallocating arrays is often better than growing arrays inside a loop,
# especially for larger simulations.


# %% Cell 6 — Element-wise version

normalized_signals_elementwise = signals ./ signals[1]

println("Normalized signals using element-wise division:")
println(normalized_signals_elementwise)

# This gives the same result as the loop above, but is shorter.


# %% Cell 7 — Log-normalized signal

log_normalized_signals = log.(normalized_signals_elementwise)

println("Log-normalized signals:")
println(log_normalized_signals)

# The dot in log.() applies log() to each element of the vector.
# This is useful for ADC fitting:
#
#     ln(S / S0) = c - ADC * b


# %% Cell 8 — Simple statistics

mean_signal = mean(signals)
std_signal = std(signals)
min_signal = minimum(signals)
max_signal = maximum(signals)

println("Signal summary:")
println("Mean = ", mean_signal)
println("Std = ", std_signal)
println("Min = ", min_signal)
println("Max = ", max_signal)


# %% Cell 9 — Dictionaries

sequence_parameters = Dict(
    "TE_ms" => 80,
    "TR_ms" => 300,
    "diffusivity" => 2.0,
    "nspins" => 2000,
)

println("Sequence parameters dictionary:")
println(sequence_parameters)

println("TE = ", sequence_parameters["TE_ms"], " ms")
println("TR = ", sequence_parameters["TR_ms"], " ms")

# A dictionary stores key-value pairs.
# It is useful for grouping parameters in a readable way.


# %% Cell 10 — Creating a small result table

header = ["bval" "signal" "normalized_signal" "log_normalized_signal"]

data = hcat(
    bvals,
    signals,
    normalized_signals_elementwise,
    log_normalized_signals,
)

table = vcat(header, data)

println("Result table:")
println(table)


# %% Cell 11 — Saving output to a file

output_dir = joinpath("docs", "data", "processed")
mkpath(output_dir)

output_file = joinpath(output_dir, "first_julia_test_signal_table.csv")

writedlm(output_file, table, ',')

println("Saved table to:")
println(output_file)

# This introduces the same output logic used later for simulation results.


# %% Cell 12 — End of script

println("Arrays, loops, statistics, and file output script completed successfully.")
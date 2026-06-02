# ============================================================
# 02_functions_and_structured_script.jl
#
# Goal:
# This script shows how to organize a Julia script using functions.
#
# This structure is used later in the simulation scripts:
#   - define constants
#   - define helper functions
#   - define a main() function
#   - call main() at the end
#
# Run from the repository root with:
#
#     julia --project=. docs/scripts/01_first_julia_tests/02_functions_and_structured_script.jl
#
# Or run cell by cell in VS Code.
# ============================================================


# %% Cell 1 — Simple function

function add_numbers(x, y)
    return x + y
end

result = add_numbers(3, 5)

println("3 + 5 = ", result)


# %% Cell 2 — Function without explicit return

function multiply_numbers(x, y)
    x * y
end

result = multiply_numbers(3, 5)

println("3 * 5 = ", result)

# In Julia, the last expression of a function is returned automatically.
# However, using return explicitly is often clearer for beginners.


# %% Cell 3 — Function with typed arguments

function compute_normalized_signal(signal::Float64, reference::Float64)
    return signal / reference
end

normalized_signal = compute_normalized_signal(0.82, 1.0)

println("Normalized signal = ", normalized_signal)

# Type annotations are optional in Julia, but they can make scientific
# scripts clearer when the expected input is known.


# %% Cell 4 — Function applied to a vector

function normalize_signal_vector(signals)
    reference = signals[1]
    return signals ./ reference
end

signals = [1.0, 0.82, 0.67, 0.55, 0.45]
normalized_signals = normalize_signal_vector(signals)

println("Signals:")
println(signals)

println("Normalized signals:")
println(normalized_signals)


# %% Cell 5 — Constants

const BVALS = [0.0, 0.5, 1.0, 1.5, 2.0]
const TE = 80
const TR = 300

println("B-values: ", BVALS)
println("TE = ", TE, " ms")
println("TR = ", TR, " ms")

# Constants are useful for simulation parameters that should not change
# during the script.


# %% Cell 6 — Helper function for displaying parameters

function print_sequence_parameters(bvals, TE, TR)
    println("Sequence parameters")
    println("-------------------")
    println("b-values: ", bvals)
    println("TE = ", TE, " ms")
    println("TR = ", TR, " ms")
end

print_sequence_parameters(BVALS, TE, TR)


# %% Cell 7 — Main function

function main()
    println("Starting structured Julia script...")

    print_sequence_parameters(BVALS, TE, TR)

    example_signals = [1.0, 0.82, 0.67, 0.55, 0.45]
    example_normalized_signals = normalize_signal_vector(example_signals)

    println("\nExample signals:")
    println(example_signals)

    println("\nExample normalized signals:")
    println(example_normalized_signals)

    println("\nStructured Julia script completed successfully.")
end


# %% Cell 8 — Script execution

main()

# In larger scripts, keeping all execution inside main() helps make the
# workflow easier to read, test, and modify.
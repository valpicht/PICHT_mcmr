# ============================================================
# 01_basic_julia_syntax.jl
#
# Goal:
# This script introduces the basic syntax of Julia.
#
# It is written with cell markers (# %%) so that it can be run
# step by step in VS Code, similarly to a Jupyter notebook.
#
# Run the full script from the repository root with:
#
#     julia --project=. docs/scripts/01_first_julia_tests/01_basic_julia_syntax.jl
#
# Or run cell by cell in VS Code.
# ============================================================


# %% Cell 1 — First print statement

println("Hello from Julia!")

# In Julia, comments start with #.
# println() prints text or values to the terminal.


# %% Cell 2 — Variables

x = 3
y = 5

println("x = ", x)
println("y = ", y)

# Julia automatically infers the type of variables.
println("Type of x: ", typeof(x))
println("Type of y: ", typeof(y))


# %% Cell 3 — Basic arithmetic

sum_xy = x + y
difference_xy = x - y
product_xy = x * y
division_xy = x / y
power_xy = x^2

println("x + y = ", sum_xy)
println("x - y = ", difference_xy)
println("x * y = ", product_xy)
println("x / y = ", division_xy)
println("x^2 = ", power_xy)


# %% Cell 4 — Strings

project_name = "MCMR neuromelanin simulation project"

println("Project name: ", project_name)
println("Type of project_name: ", typeof(project_name))

# String interpolation is often cleaner than concatenation.
println("This project is called: $project_name")


# %% Cell 5 — Booleans and comparisons

is_ready = true

println("Is the script ready? ", is_ready)
println("Type of is_ready: ", typeof(is_ready))

println("Is x smaller than y? ", x < y)
println("Is x equal to y? ", x == y)
println("Is x different from y? ", x != y)


# %% Cell 6 — If / else statements

if x < y
    println("x is smaller than y.")
elseif x == y
    println("x is equal to y.")
else
    println("x is larger than y.")
end


# %% Cell 7 — Vectors

bvals = [0.0, 0.5, 1.0, 1.5, 2.0]

println("b-values:")
println(bvals)

println("First b-value: ", bvals[1])
println("Last b-value: ", bvals[end])
println("Number of b-values: ", length(bvals))

# Important:
# Julia uses 1-based indexing, not 0-based indexing.
# The first element is bvals[1], not bvals[0].


# %% Cell 8 — Element-wise operations

signals = [1.0, 0.82, 0.67, 0.55, 0.45]

normalized_signals = signals ./ signals[1]

println("Original signals:")
println(signals)

println("Normalized signals:")
println(normalized_signals)

# The dot . means element-wise operation.
# signals ./ signals[1] divides every element of signals by signals[1].


# %% Cell 9 — End of script

println("Basic Julia syntax script completed successfully.")
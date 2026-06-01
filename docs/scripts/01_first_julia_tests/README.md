# 01 — First Julia tests

This folder contains introductory Julia scripts used before running the MCMR simulations.

The goal is to introduce the basic Julia concepts needed later in the project, without using `MCMRSimulator` yet.

The scripts are written as `.jl` files with cell markers (`# %%`), so they can be run step by step in VS Code like lightweight notebooks.

## Workflow

The workflow is:

1. learn basic Julia syntax;
2. organize code into functions;
3. use a `main()` function;
4. work with arrays and loops;
5. compute simple statistics;
6. save a small output table.

These concepts are reused later for MCMR simulation scripts, especially when looping over b-values, echo times, substrates, and simulation outputs.

## How to run

From the repository root, run for example:

```bash
julia +1.11 --project=. docs/scripts/01_first_julia_tests/01_basic_julia_syntax.jl
```

The scripts can also be run cell by cell in VS Code.

## Scripts

### `01_basic_julia_syntax.jl`

Introduces the most basic Julia syntax.

This script covers:

- printing with `println`;
- comments;
- variables;
- basic arithmetic;
- strings;
- string interpolation;
- booleans;
- comparisons;
- `if` / `elseif` / `else` statements;
- vectors;
- 1-based indexing;
- element-wise operations.

Run with:

```bash
julia +1.11 --project=. docs/scripts/01_first_julia_tests/01_basic_julia_syntax.jl
```

### `02_functions_and_structured_script.jl`

Introduces how to organize a Julia script using functions.

This script covers:

- defining simple functions;
- explicit `return`;
- implicit return values;
- typed function arguments;
- applying a function to a vector;
- using constants for simulation parameters;
- defining helper functions;
- organizing execution inside a `main()` function.

This structure is used later in the simulation scripts to keep code readable and reproducible.

Run with:

```bash
julia +1.11 --project=. docs/scripts/01_first_julia_tests/02_functions_and_structured_script.jl
```

### `03_arrays_loops_and_statistics.jl`

Introduces arrays, loops, simple statistics, dictionaries, and file output.

This script covers:

- vectors;
- loops over one vector;
- loops over two vectors using `zip`;
- preallocating arrays;
- element-wise operations;
- log-normalized signals;
- basic statistics with `Statistics`;
- dictionaries for parameter storage;
- creating a small result table;
- saving a CSV-like output file with `DelimitedFiles`.

Run with:

```bash
julia +1.11 --project=. docs/scripts/01_first_julia_tests/03_arrays_loops_and_statistics.jl
```

Output:

```text
docs/data/processed/first_julia_test_signal_table.csv
```

## Notes

These scripts do not use `MCMRSimulator` yet.

The next folder introduces the first MCMR examples:

```text
docs/scripts/02_first_mcmr_tests/
```

Those scripts move from basic Julia syntax to actual MRI simulation objects such as geometries, sequences, simulations, readouts, and snapshots.
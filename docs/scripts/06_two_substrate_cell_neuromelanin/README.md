# 06 — Two-substrate cell + neuromelanin simulation

This folder contains scripts for simulations with two substrates:

1. a cell substrate;
2. a neuromelanin substrate.

Both substrates are loaded from SWC files and converted into overlapping-sphere geometries. The goal is to reproduce the main project setup where local relaxation properties are assigned to different compartments and gradient-echo signals are extracted separately.

This workflow is more specific than the previous single-substrate DWI scripts. It is designed for the project’s main neuromelanin-sensitive MRI simulations.

## Required input

The scripts expect two SWC files in:

```text
substrate_generation/example_outputs/
```

The default files are:

```text
substrate_generation/example_outputs/cell_soma_r19_one_branche_neurons_list_new.swc
substrate_generation/example_outputs/neuromelanin_sphere_r15_neurons_list_new.swc
```

If the files are missing, copy them from the substrate generation output, for example:

```bash
cp ~/project/Grow_neurons/results/cell_soma_r19_one_branche_neurons_list_new.swc \
   substrate_generation/example_outputs/

cp ~/project/Grow_neurons/results/neuromelanin_sphere_r15_neurons_list_new.swc \
   substrate_generation/example_outputs/
```

In the scripts, the input files are selected with:

```julia
const CELL_SWC_NAME = "cell_soma_r19_one_branche_neurons_list_new"
const NM_SWC_NAME = "neuromelanin_sphere_r15_neurons_list_new"
```

The `.swc` extension is optional.

## Geometry convention

The combined geometry is always built in this order:

```julia
geometry = [
    cell_geometry,
    neuromelanin_geometry,
]
```

Therefore:

```text
geometry[1] = cell
geometry[2] = neuromelanin
```

This is important because compartment-specific readouts use `geometry_index`:

```julia
Subset(inside = true, geometry_index = 1)  # inside cell
Subset(inside = false, geometry_index = 1) # outside cell
Subset(inside = true, geometry_index = 2)  # inside neuromelanin
```

## Relaxation convention

The scripts use time in milliseconds. Therefore, relaxation rates should be given in `ms⁻¹`.

Examples:

```text
T1 = 1000 ms  -> R1 = 1 / 1000 ms⁻¹ = 0.001 ms⁻¹
T2 = 3 ms     -> R2 = 1 / 3 ms⁻¹ ≈ 0.333 ms⁻¹
```

The simulations use:

```julia
global_R1
global_R2
cell_R1_inside
nm_R1_inside
```

The local R1 inside each geometry is treated as an additional local contribution:

```text
local R1 = global_R1
         + cell_R1_inside if inside cell
         + nm_R1_inside if inside neuromelanin
```

## Utility file

### `utils_cell_neuromelanin.jl`

This file contains reusable helper functions for the two-substrate workflow.

It includes functions to:

- build relative paths to SWC files;
- read SWC files;
- extract sphere positions and radii;
- create the cell geometry;
- create the neuromelanin geometry;
- assign local `R1_inside` values;
- create a common bounding box around both substrates;
- generate random spins;
- count spins inside/outside each compartment;
- check expected local R1 values at selected positions;
- create a `GradientEcho` sequence;
- run gradient-echo simulations;
- extract total and compartment-wise signals;
- save CSV outputs.

This avoids repeating the same code in every script and avoids absolute paths such as:

```text
/home/valentine/project/...
```

## Scripts

### `01_check_two_substrate_geometry_and_R1.jl`

Checks that the two-substrate setup is correct before running the full gradient-echo simulation.

The script:

1. loads the cell SWC substrate;
2. loads the neuromelanin SWC substrate;
3. creates two overlapping-sphere geometries;
4. assigns local R1 contributions;
5. creates a common bounding box;
6. generates random spins;
7. counts spins inside and outside each geometry;
8. checks expected local R1 values at selected positions;
9. saves a local R1 check CSV file.

Run with:

```bash
julia +1.11 --project=. docs/scripts/06_two_substrate_cell_neuromelanin/01_check_two_substrate_geometry_and_R1.jl
```

Output:

```text
docs/data/processed/cell_r19_nm_r15_local_R1_check.csv
```

The output CSV contains:

```text
label
x, y, z
inside_cell
inside_neuromelanin
global_R1
cell_R1_inside
nm_R1_inside
expected_local_R1
```

This script is mainly a validation step. It should be run before the full gradient-echo simulation.

### `02_run_cell_neuromelanin_gradient_echo.jl`

Runs a gradient-echo simulation with the cell and neuromelanin substrates.

The script:

1. loads the cell substrate;
2. loads the neuromelanin substrate;
3. creates the combined geometry;
4. assigns local R1 values;
5. sets global R2;
6. generates random spins in a common bounding box;
7. runs a `GradientEcho` simulation;
8. extracts compartment-wise signals;
9. saves the results to CSV.

Run with:

```bash
julia +1.11 --project=. docs/scripts/06_two_substrate_cell_neuromelanin/02_run_cell_neuromelanin_gradient_echo.jl
```

Output:

```text
docs/data/processed/cell_r19_nm_r15_gradient_echo_compartment_signals.csv
```

The output CSV contains one row per neuromelanin R1 condition.

Main user parameters are defined at the top of the script:

```julia
const CELL_SWC_NAME = "cell_soma_r19_one_branche_neurons_list_new"
const NM_SWC_NAME = "neuromelanin_sphere_r15_neurons_list_new"

const OUTPUT_PREFIX = "cell_r19_nm_r15"

const NSPINS = 50_000
const BBOX_MARGIN = 2.0
const RANDOM_SEED = 1234

const TE = 3
const TR = 25

const DIFFUSIVITY = 1.0
const SKIP_TR = 20

const GLOBAL_R1 = 0.0
const CELL_R1_INSIDE = 0.0005
const NM_R1_INSIDE_VALUES = [
    0.0005,
    0.0010,
    1.0000,
]

const GLOBAL_R2 = 1.0 / 3.0

const CELL_PERMEABILITY = 0.0
const NM_PERMEABILITY = Inf
```

The extracted signals are:

```text
total_signal
inside_cell_signal
outside_cell_signal
inside_neuromelanin_signal
```

The corresponding spin counts are also saved:

```text
total_nspins
inside_cell_nspins
outside_cell_nspins
inside_neuromelanin_nspins
```

### `03_plot_cell_neuromelanin_signals.jl`

Plots the compartment-wise gradient-echo signals generated by the simulation script.

Run with:

```bash
julia +1.11 --project=. docs/scripts/06_two_substrate_cell_neuromelanin/03_plot_cell_neuromelanin_signals.jl
```

Input:

```text
docs/data/processed/cell_r19_nm_r15_gradient_echo_compartment_signals.csv
```

Outputs:

```text
docs/figures/cell_r19_nm_r15_gradient_echo_absolute_compartment_signals.png
docs/figures/cell_r19_nm_r15_gradient_echo_normalized_compartment_signals.png
docs/figures/cell_r19_nm_r15_gradient_echo_neuromelanin_signal.png
```

The plots show:

1. absolute compartment-wise signals;
2. signals normalized to the first neuromelanin R1 condition;
3. neuromelanin-only signal, both absolute and normalized.

## Interpretation guide

The total signal is the combined gradient-echo signal from all spins.

The inside cell signal corresponds to spins inside the cell geometry:

```julia
Subset(inside = true, geometry_index = 1)
```

The outside cell signal corresponds to spins outside the cell geometry:

```julia
Subset(inside = false, geometry_index = 1)
```

The inside neuromelanin signal corresponds to spins inside the neuromelanin geometry:

```julia
Subset(inside = true, geometry_index = 2)
```

If changing `nm_R1_inside` mainly affects the inside-neuromelanin signal, then the neuromelanin compartment is driving the observed signal change.

If the total signal also changes, this suggests that local neuromelanin relaxation contributes measurably to the combined gradient-echo signal.

## Notes

This folder focuses on the two-substrate cell + neuromelanin workflow.

It is intended to be the bridge between the earlier validation scripts and the final project simulations.

The current scripts are designed to be reusable: to test another cell or neuromelanin substrate, copy the corresponding SWC files into:

```text
substrate_generation/example_outputs/
```

and update:

```julia
const CELL_SWC_NAME = "..."
const NM_SWC_NAME = "..."
```

in the simulation scripts.
# 03 — Triangular mesh tests

This folder contains the first simulations using triangular mesh substrates.

The goal is to move from the simple cylinder geometries used in the MCMR tutorial to custom mesh-based substrates generated during the project.

The scripts are written with cell markers (`# %%`), so they can be run step by step in VS Code like lightweight notebooks.

## Required input

The scripts expect at least one example triangular mesh file in:

```text
substrate_generation/example_outputs/
```

By default, the scripts use:

```text
substrate_generation/example_outputs/cell_sphere_r19_neurons_list_new_mesh.ply
```

If the file is missing, copy it from the substrate generation output, for example:

```bash
cp ~/project/Grow_neurons/results/cell_sphere_r19_neurons_list_new_mesh.ply \
   substrate_generation/example_outputs/
```

Other `.ply` files can also be used. In the scripts, the mesh is selected with:

```julia
const MESH_NAME = "cell_sphere_r19_neurons_list_new_mesh"
```

The `.ply` extension is optional. For example, both of these are valid:

```julia
const MESH_NAME = "cell_sphere_r19_neurons_list_new_mesh"
const MESH_NAME = "cell_sphere_r19_neurons_list_new_mesh.ply"
```

The mesh-loading logic is defined in:

```text
utils_mesh.jl
```

## Utility file

### `utils_mesh.jl`

This file contains helper functions reused by the mesh scripts.

It includes functions to:

- build relative paths to mesh files;
- load `.ply` meshes;
- print mesh coordinate ranges;
- create a bounding box around the mesh;
- create a basic DWI sequence.

This avoids repeating the same code in every script and avoids absolute paths such as:

```text
/home/valentine/project/...
```

## Scripts

### `01_load_and_plot_mesh.jl`

Loads a `.ply` triangular mesh, prints basic information about its vertices, and plots the geometry.

Run with:

```bash
julia +1.11 --project=. docs/scripts/03_triangular_mesh_tests/01_load_and_plot_mesh.jl
```

This script is mainly a smoke test to check that the mesh can be loaded and inspected correctly.

### `02_mesh_bounding_box_and_subsets.jl`

Creates a bounding box around the mesh, generates random spins, and counts how many spins are inside the mesh.

Run with:

```bash
julia +1.11 --project=. docs/scripts/03_triangular_mesh_tests/02_mesh_bounding_box_and_subsets.jl
```

This script is useful to check that the mesh can be used to separate inside and outside compartments.

### `03_run_simple_dwi_on_mesh.jl`

Runs a first diffusion-weighted simulation on the triangular mesh and separates the signal into inside and outside contributions.

Run with:

```bash
julia +1.11 --project=. docs/scripts/03_triangular_mesh_tests/03_run_simple_dwi_on_mesh.jl
```

This script is the mesh equivalent of the simple DWI simulation used earlier with the cylinder geometry.

### `04_signal_vs_b_on_mesh.jl`

Runs simulations over several b-values and saves normalized inside/outside signal curves.

Run with:

```bash
julia +1.11 --project=. docs/scripts/03_triangular_mesh_tests/04_signal_vs_b_on_mesh.jl
```

Outputs:

```text
docs/data/processed/mesh_signal_vs_b_data.csv
docs/figures/mesh_inside_signal_vs_b.png
docs/figures/mesh_outside_signal_vs_b.png
```

### `05_signal_vs_b_TE_and_ADC.jl`

Runs simulations over several b-values and echo times, computes `ln(S/S0)`, fits a simple ADC model, and saves summary plots.

Run with:

```bash
julia +1.11 --project=. docs/scripts/03_triangular_mesh_tests/05_signal_vs_b_TE_and_ADC.jl
```

Outputs:

```text
docs/data/processed/mesh_ln_signal_vs_b_TE_data.csv
docs/data/processed/mesh_adc_vs_TE_data.csv
docs/figures/mesh_inside_ln_signal_vs_b_TE_fit.png
docs/figures/mesh_outside_ln_signal_vs_b_TE_fit.png
docs/figures/mesh_adc_vs_TE.png
```

## ADC model

The ADC fit uses the simple log-linear model:

```text
ln(S/S0) = c - ADC * b
```

where:

- `S` is the signal at a given b-value;
- `S0` is the signal at `b = 0`;
- `c` is the fitted intercept;
- `ADC` is the apparent diffusion coefficient.

This is used here as a first validation step to check whether the mesh-based simulations produce reasonable diffusion-weighted signal decay.

## Notes

These scripts are still early validation scripts. They are not the final neuromelanin simulations.

The next steps in the repository introduce:

- overlapping spheres;
- comparison between mesh and overlapping-sphere substrates;
- SWC-based substrates;
- local R1/T1 effects;
- global R2 effects;
- neuromelanin-containing substrates.
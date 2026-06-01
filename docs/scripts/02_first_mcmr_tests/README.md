# 02 — First MCMR tests

This folder contains the first scripts using `MCMRSimulator`.

The goal is to reproduce the main logic of the official MCMR Julia tutorial in a cleaner, step-by-step format. These examples are written as `.jl` scripts with cell markers (`# %%`), so they can be run cell by cell in VS Code like lightweight notebooks.

The workflow is:

1. create and plot a simple geometry;
2. define a diffusion-weighted MRI sequence;
3. create a `Simulation` object;
4. run simple signal readouts;
5. split the signal into inside/outside compartments;
6. visualize signal evolution and spin snapshots.

## Relation to the official MCMR tutorial

These scripts are based on the first examples from the MCMR Julia tutorial, but reorganized into smaller files.

The main difference is that each script focuses on one concept only, so the workflow is easier to follow and easier to debug.

## Scripts

### `01_create_and_plot_geometry.jl`

Creates a simple repeated-cylinder geometry and plots it.

This introduces:

- loading `MCMRSimulator`;
- creating a `Cylinders` geometry;
- using `PlotPlane`;
- plotting a geometry with CairoMakie.

Run with:

```bash
julia +1.11 --project=. docs/scripts/02_first_mcmr_tests/01_create_and_plot_geometry.jl
```

### `02_create_sequence_and_simulation.jl`

Creates a diffusion-weighted MRI sequence using `MRIBuilder`, then combines the sequence and geometry into a `Simulation` object.

This introduces:

- loading `MRIBuilder`;
- creating a `DWI` sequence;
- plotting the sequence with `plot_sequence`;
- defining basic simulation parameters such as `R1`, `R2`, diffusivity, off-resonance, and geometry.

Run with:

```bash
julia +1.11 --project=. docs/scripts/02_first_mcmr_tests/02_create_sequence_and_simulation.jl
```

### `03_simple_readouts_and_subsets.jl`

Runs simple signal readouts and separates the signal into inside-cylinder and outside-cylinder contributions.

This introduces:

- `readout(nspins, simulation)`;
- `skip_TR`;
- `Subset(inside = true)`;
- `Subset(inside = false)`;
- extracting longitudinal and transverse signal components.

Run with:

```bash
julia +1.11 --project=. docs/scripts/02_first_mcmr_tests/03_simple_readouts_and_subsets.jl
```

### `04_signal_evolution_and_snapshots.jl`

Visualizes signal evolution over repeated TRs and returns spin snapshots.

This introduces:

- readouts over several TRs with `nTR`;
- signal evolution over custom time points;
- trajectory snapshots with `return_snapshot = true`;
- plotting snapshots in different ways.

Run with:

```bash
julia +1.11 --project=. docs/scripts/02_first_mcmr_tests/04_signal_evolution_and_snapshots.jl
```

## Notes

These scripts use the simple `Cylinders` geometry from the MCMR tutorial.

The next folder introduces triangular mesh substrates generated during the project:

```text
docs/scripts/03_triangular_mesh_tests/
```

The first mesh scripts move from simple built-in geometries to project-specific `.ply` triangular meshes.
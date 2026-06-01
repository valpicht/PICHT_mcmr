# 04 — Overlapping spheres

This folder contains the first scripts using the overlapping-sphere substrate representation implemented during the project.

The goal is to move from triangular mesh substrates to a sphere-based representation where multiple spheres can overlap and together define a single biological-like object.

This is an important step because SWC files naturally describe neurite-like geometries as connected points with radii. These points can be interpreted as spheres, and allowing overlap makes it possible to approximate a continuous cellular substrate.

## Required input

The scripts expect at least one SWC file in:

```text
substrate_generation/example_outputs/
```

For example:

```text
substrate_generation/example_outputs/one_soma_list_new.swc
```

If the file is missing, copy it from the substrate generation output, for example:

```bash
cp ~/project/Grow_neurons/results/one_soma_list_new.swc \
   substrate_generation/example_outputs/
```

The directional validation script also expects a long cylinder-like SWC file, for example:

```text
substrate_generation/example_outputs/cylinder_soma_neurons_list_new.swc
```

If the file is missing, copy it with:

```bash
cp ~/project/Grow_neurons/results/cylinder_soma_neurons_list_new.swc \
   substrate_generation/example_outputs/
```

Other `.swc` files can also be used. In the scripts, the substrate is selected with:

```julia
const SWC_NAME = "one_soma_list_new"
```

or, for the long-cylinder validation:

```julia
const SWC_NAME = "cylinder_soma_neurons_list_new"
```

The `.swc` extension is optional.

## Utility file

### `utils_overlapping.jl`

This file contains helper functions reused by the overlapping-sphere scripts.

It includes functions to:

- build relative paths to SWC files;
- read SWC files;
- extract sphere positions and radii;
- create MCMR `Spheres` geometries;
- create an MCMR `Spheres` geometry with `overlapping = true`;
- create a bounding box around the spheres;
- generate random spins;
- count inside/outside spins;
- create DWI sequences and simulations;
- run b-value sweeps;
- normalize signals by `S(0)`;
- compute `ln(S/S0)`;
- save signal results to CSV.

This avoids repeating the same code in every script and avoids absolute paths such as:

```text
/home/valentine/project/...
```

## Scripts

### `01_create_overlapping_spheres.jl`

Loads an SWC file, extracts sphere positions and radii, creates an overlapping-sphere geometry, and checks how many randomly generated spins are inside or outside the geometry.

Run with:

```bash
julia +1.11 --project=. docs/scripts/04_overlapping_spheres/01_create_overlapping_spheres.jl
```

This script is mainly a smoke test to verify that the SWC-to-overlapping-spheres conversion works.

### `02_run_dwi_overlapping_spheres.jl`

Runs a first diffusion-weighted simulation on the overlapping-sphere substrate.

Run with:

```bash
julia +1.11 --project=. docs/scripts/04_overlapping_spheres/02_run_dwi_overlapping_spheres.jl
```

Outputs:

```text
docs/data/processed/overlapping_spheres_signal_vs_b_data.csv
docs/figures/overlapping_spheres_inside_signal_vs_b.png
docs/figures/overlapping_spheres_outside_signal_vs_b.png
```

### `03_compare_non_overlapping_vs_overlapping.jl`

Compares the same SWC-derived sphere substrate with:

```julia
overlapping = false
```

and:

```julia
overlapping = true
```

This script uses the same sphere positions, radii, spin distribution, b-values, and DWI sequence parameters for both cases. It is intended to validate the effect of the `overlapping` option itself.

Run with:

```bash
julia +1.11 --project=. docs/scripts/04_overlapping_spheres/03_compare_non_overlapping_vs_overlapping.jl
```

Outputs:

```text
docs/data/processed/non_overlapping_vs_overlapping_signal_vs_b_data.csv
docs/figures/non_overlapping_vs_overlapping_inside_signal_vs_b.png
docs/figures/non_overlapping_vs_overlapping_outside_signal_vs_b.png
```

### `04_long_cylinder_free_vs_restricted_diffusion.jl`

Validates the directional diffusion behavior of an overlapping-sphere substrate using a long cylinder-like SWC geometry.

The cylinder is aligned along the x-axis. The expected qualitative behavior is:

- diffusion along x is less restricted;
- diffusion along y and z is more restricted;
- diffusion weighting along x should produce stronger signal attenuation;
- the apparent ADC should be higher along x than along y or z.

Run with:

```bash
julia +1.11 --project=. docs/scripts/04_overlapping_spheres/04_long_cylinder_free_vs_restricted_diffusion.jl
```

Outputs:

```text
docs/data/processed/long_cylinder_directional_signal.csv
docs/data/processed/long_cylinder_directional_ADC_summary.csv
docs/figures/long_cylinder_directional_signal_vs_b.png
docs/figures/long_cylinder_directional_ln_signal_fit.png
docs/figures/long_cylinder_directional_ADC_summary.png
```

This script provides a qualitative validation of the overlapping-sphere implementation. It checks that a long cylinder-like substrate shows higher apparent diffusion along its long axis than along transverse directions.

## Important note about single-sphere substrates

The current example file:

```text
one_soma_list_new.swc
```

contains only one sphere. Therefore, `overlapping = false` and `overlapping = true` are expected to give very similar or identical results for this specific example.

The difference between the two options becomes meaningful for substrates containing multiple spheres that geometrically overlap.

## Interpretation of the long-cylinder validation

In diffusion MRI, stronger diffusion along the direction of the diffusion gradient produces stronger signal attenuation.

Therefore, for a cylinder aligned along x:

```text
free or less restricted diffusion along x
→ stronger signal attenuation in x
→ steeper negative slope in ln(S/S0)
→ higher ADC along x
```

The expected qualitative result is:

```text
ADC_x > ADC_y and ADC_x > ADC_z
```

This validates that the overlapping-sphere geometry can reproduce direction-dependent diffusion behavior.

## Notes

This folder focuses on the overlapping-sphere representation and validation of the `overlapping = true` option.

The next step in the project is to move from validation tests to more general simulation scripts that can run on arbitrary SWC substrates, and then to combined cell + neuromelanin simulations.
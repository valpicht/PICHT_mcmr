# 07 — Cell + neuromelanin R1 sweep

This folder contains the scripts used to run a parameter sweep on the local relaxation rate assigned to the neuromelanin compartment.

The goal is to simulate how the neuromelanin-sensitive MRI contrast ratio changes when the total neuromelanin R1 value is varied while keeping the rest of the simulation fixed.

This step follows the supervisor’s request to test several neuromelanin R1 values, including values around the estimated value of approximately:

```text
R1_NM ≈ 50 s^-1
```

The contrast ratio is computed as:

```text
CR = (S_i - S_0) / S_0
```

where:

- `S_i` is the signal for a given simulated neuromelanin R1 value;
- `S_0` is the baseline signal where the neuromelanin R1 is equal to the background R1.

In this workflow, the baseline is:

```text
R1_NM = R1_background = 1 s^-1
```

## Scientific purpose

The simulation tests whether changing the local R1 assigned to the neuromelanin compartment can modify the simulated MRI signal and the resulting contrast ratio.

The biological motivation comes from neuromelanin-sensitive MRI, where contrast is often quantified using a signal ratio between a neuromelanin-containing region and a reference region. In this simplified simulation, the contrast ratio is computed relative to a baseline simulation where neuromelanin has the same R1 as the background.

This is not yet a full biological model of neuromelanin MRI contrast. It is a controlled simulation experiment designed to test the effect of one parameter:

```text
R1_NM
```

while keeping the geometry, sequence parameters, global R1, global R2, diffusivity, and spin distribution fixed.

## Required input

The scripts expect two SWC files in:

```text
substrate_generation/example_outputs/
```

One file represents the cell substrate, and one file represents the neuromelanin substrate.

Example expected files:

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

The names used by the simulation are defined at the top of:

```text
01_run_cell_nm_R1_sweep.jl
```

with:

```julia
const CELL_SWC_NAME = "cell_soma_r19_one_branche_neurons_list_new"
const NM_SWC_NAME = "neuromelanin_sphere_r15_neurons_list_new"
```

The `.swc` extension is optional.

## Utility file

### `utils_cell_nm_R1_sweep.jl`

This file contains reusable helper functions for the cell + neuromelanin R1 sweep.

It includes functions to:

- convert relaxation rates between `s^-1` and `ms^-1`;
- build relative paths to SWC files;
- load SWC files as sphere positions and radii;
- create the cell geometry;
- create the neuromelanin geometry;
- assign additional local R1 to the neuromelanin compartment;
- create a combined bounding box around both substrates;
- generate a reproducible spin distribution;
- count spins inside the cell and neuromelanin compartments;
- create a gradient-echo sequence;
- run compartment-wise readouts;
- compute contrast ratios;
- save the R1 sweep results to CSV;
- save a local R1 check table to CSV.

The utility file avoids hard-coded absolute paths such as:

```text
/home/valentine/project/...
```

and makes the workflow easier to reuse on another computer.

## Scripts

### `01_run_cell_nm_R1_sweep.jl`

Runs the cell + neuromelanin gradient-echo simulations while sweeping the total R1 value assigned to the neuromelanin compartment.

Run with:

```bash
julia +1.11 --project=. docs/scripts/07_cell_neuromelanin_R1_sweep/01_run_cell_nm_R1_sweep.jl
```

The script:

1. loads the cell SWC substrate;
2. loads the neuromelanin SWC substrate;
3. creates a two-substrate geometry:
   - geometry index 1: cell;
   - geometry index 2: neuromelanin;
4. applies a global/background R1 to the whole simulation;
5. applies an additional local R1 contribution inside the neuromelanin geometry;
6. applies a global R2;
7. runs a gradient-echo readout;
8. extracts compartment-wise signals;
9. computes contrast ratios relative to the baseline;
10. saves the results progressively.

Main parameters are defined at the top of the script:

```julia
const CELL_SWC_NAME = "cell_soma_r19_one_branche_neurons_list_new"
const NM_SWC_NAME = "neuromelanin_sphere_r15_neurons_list_new"

const OUTPUT_PREFIX = "cell_nm_R1_sweep"

const R1_BACKGROUND_S = 1.0
const R2_GLOBAL_S = 1.0 / 3.0

const R1_NM_TOTAL_VALUES_S = default_R1_NM_values_s()

const TE = 3
const TR = 25

const DIFFUSIVITY = 1.0
const NSPINS = 50_000
const RANDOM_SEED = 1234
const BBOX_MARGIN = 2.0
const SKIP_TR = 20
const MANUAL_TIMESTEP = nothing
```

The default R1 sweep is:

```text
1, 10, 20, 30, ..., 100, 200, 300, ..., 1000 s^-1
```

This is defined in:

```julia
default_R1_NM_values_s()
```

If the full sweep takes too long, a coarser sweep can be used instead:

```julia
const R1_NM_TOTAL_VALUES_S = coarse_R1_NM_values_s()
```

The coarse sweep is:

```text
1, 20, 40, 60, 80, 100, 200, 400, 600, 800, 1000 s^-1
```

#### Outputs

The script creates:

```text
docs/data/processed/cell_nm_R1_sweep_R1_sweep.csv
docs/data/processed/cell_nm_R1_sweep_local_R1_check.csv
```

The main sweep CSV contains:

```text
R1_background_s
R1_NM_total_s
R1_NM_additional_s
R2_global_s
TE
TR
n_total
n_inside_cell
n_outside_cell
n_inside_nm
n_outside_nm
total_signal
inside_cell_signal
outside_cell_signal
inside_nm_signal
outside_nm_signal
CR_total
CR_inside_cell
CR_outside_cell
CR_inside_nm
CR_outside_nm
```

The local R1 check CSV contains simple diagnostic rows indicating whether selected points are inside the cell and/or neuromelanin geometry and what R1 value is expected there.

### `02_plot_cell_nm_R1_sweep_CR.jl`

Plots the results generated by:

```text
01_run_cell_nm_R1_sweep.jl
```

Run with:

```bash
julia +1.11 --project=. docs/scripts/07_cell_neuromelanin_R1_sweep/02_plot_cell_nm_R1_sweep_CR.jl
```

The script reads:

```text
docs/data/processed/cell_nm_R1_sweep_R1_sweep.csv
```

and creates the following figures:

```text
docs/figures/cell_nm_R1_sweep_CR_total_vs_R1_NM.png
docs/figures/cell_nm_R1_sweep_CR_compartments_vs_R1_NM.png
docs/figures/cell_nm_R1_sweep_raw_signals_vs_R1_NM.png
docs/figures/cell_nm_R1_sweep_CR_total_zoom_0_to_120_s.png
```

#### Figure 1 — Total contrast ratio

```text
cell_nm_R1_sweep_CR_total_vs_R1_NM.png
```

Shows the total simulated contrast ratio as a function of total neuromelanin R1.

This is the main figure that directly answers the requested analysis:

```text
x-axis: R1_NM total (s^-1)
y-axis: CR = (S_i - S_0) / S_0
```

A vertical dashed line marks the estimated neuromelanin R1 value:

```text
R1_NM ≈ 50 s^-1
```

#### Figure 2 — Compartment-wise contrast ratio

```text
cell_nm_R1_sweep_CR_compartments_vs_R1_NM.png
```

Shows the contrast ratio separately for:

- total signal;
- inside cell signal;
- outside cell signal;
- inside neuromelanin signal.

This figure is useful for interpreting which compartment drives the total contrast change.

#### Figure 3 — Raw signals

```text
cell_nm_R1_sweep_raw_signals_vs_R1_NM.png
```

Shows the raw gradient-echo signal magnitudes as a function of total neuromelanin R1.

This figure is mainly a diagnostic plot. It helps check that the contrast ratio trends are not caused by unexpected numerical artefacts.

#### Figure 4 — Zoom near estimated R1_NM

```text
cell_nm_R1_sweep_CR_total_zoom_0_to_120_s.png
```

Shows the total contrast ratio in the lower R1 range:

```text
0–120 s^-1
```

This plot is useful because the estimated neuromelanin R1 is around:

```text
50 s^-1
```

and therefore this range is the most relevant for interpretation.

## Contrast ratio definition

The contrast ratio is computed as:

```text
CR = (S_i - S_0) / S_0
```

where:

- `S_i` is the signal for a given neuromelanin R1 value;
- `S_0` is the baseline signal;
- the baseline is the simulation where `R1_NM = R1_background`.

In the default configuration:

```text
R1_background = 1 s^-1
R1_NM baseline = 1 s^-1
```

so:

```text
S_0 = signal at R1_NM = 1 s^-1
```

The same formula is applied to total and compartment-wise signals.

## Important note about R1 assignment

The global/background R1 is applied to the whole simulation.

The neuromelanin geometry receives an additional local R1 contribution so that the total R1 inside neuromelanin equals the target value:

```text
R1_NM_total = R1_background + R1_NM_additional
```

Therefore:

```text
R1_NM_additional = R1_NM_total - R1_background
```

For example, if:

```text
R1_background = 1 s^-1
R1_NM_total = 50 s^-1
```

then:

```text
R1_NM_additional = 49 s^-1
```

This additional value is converted to `ms^-1` before being passed to the MCMR geometry.

## Important note about units

The user-facing parameters in the script are written in `s^-1`:

```julia
const R1_BACKGROUND_S = 1.0
const R2_GLOBAL_S = 1.0 / 3.0
const R1_NM_TOTAL_VALUES_S = default_R1_NM_values_s()
```

Inside the utility functions, these rates are converted to `ms^-1` because TE and TR are expressed in milliseconds.

The conversion is:

```text
rate_ms = rate_s / 1000
```

Be careful with `R2_GLOBAL_S`.

For example:

```text
1/3 s^-1  = 0.000333 ms^-1
1/3 ms^-1 = 333.333 s^-1
```

These are very different values.

If the intended global R2 is:

```text
R2 = 1/3 ms^-1
```

then the script should use:

```julia
const R2_GLOBAL_S = 1000.0 / 3.0
```

If the intended global R2 is:

```text
R2 = 1/3 s^-1
```

then the current default is correct:

```julia
const R2_GLOBAL_S = 1.0 / 3.0
```

This should be checked before interpreting the results quantitatively.

## Example interpretation

In the current sweep, the total contrast ratio changes with the assigned neuromelanin R1.

The zoomed plot around the estimated value is useful for presentation because it shows the simulated contrast near:

```text
R1_NM ≈ 50 s^-1
```

The compartment-wise plot helps interpret the origin of the contrast change. In the current output, the outside-cell signal changes very little, while the inside-cell and inside-neuromelanin compartments show stronger variation with R1_NM.

This suggests that the framework can be used to decompose the simulated NM-MRI signal into compartment-specific contributions.

## Recommended figures for presentation

For the main presentation, use:

```text
docs/figures/cell_nm_R1_sweep_CR_total_zoom_0_to_120_s.png
```

This is the clearest figure for showing the effect around the estimated neuromelanin R1 value.

As a second or backup figure, use:

```text
docs/figures/cell_nm_R1_sweep_CR_compartments_vs_R1_NM.png
```

This figure explains which compartment contributes to the total contrast trend.

Keep the raw signal plot as a diagnostic or backup figure:

```text
docs/figures/cell_nm_R1_sweep_raw_signals_vs_R1_NM.png
```

## Notes and limitations

This is a simplified simulation.

Important limitations:

- the geometry is simplified;
- the spin count may need to be increased for more stable results;
- repeats should be added to estimate Monte Carlo variability;
- the chosen global R2 value must be checked carefully;
- the model does not yet include magnetization transfer effects;
- the result should be interpreted as a simulation trend, not as a final quantitative prediction of in vivo NM-MRI contrast.

Recommended next checks:

1. repeat the sweep with several random seeds;
2. verify the intended R2 unit;
3. compare results for different cell/NM geometries;
4. test whether the trend is robust to spin count;
5. later add magnetization transfer effects.
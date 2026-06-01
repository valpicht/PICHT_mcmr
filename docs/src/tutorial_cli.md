# [Tutorial](@id tutorial_cli)
This tutorial will walk through an example of modelling the MRI signal evolution for a diffusion-weighted sequence.
The spins in this simulation will be constrained by regularly packed cylinders.
This tutorial will use the command line interface, which we assume is available through the `mcmr` command (see [installation instructions](@ref installation)).
If you would prefer to interact with MCMRSimulator in Julia, you can find a tutorial doing the same analysis [here](@ref tutorial_julia).

The command line tool requires one or more MRI sequence definitions as input.
These can be supplied as [pulseq](http://pulseq.github.io/) files, such as those produced by [MRIBuilder.jl](https://open.win.ox.ac.uk/pages/ndcn0236/mribuilder.jl/dev/).

In general, running a simulation will consist of the following two steps:
- Creating a geometry using one or more calls to `mcmr geometry` ([full description](@ref geometry)).
- Running the actual simulation using `mcmr run`.
We will look through each of these steps below.

## Defining the geometry
First we will define a geometry formed of regularly packed axons.
This is represented by a single cylinder pointing in the z-direction with a radius of 1 micrometer that repeats itself every 2.5 micrometer (in both the x-, and y-direction).
```bash
mcmr geometry create cylinders 1 geometry.json --radius 1 --repeats 2.5,2.5
```
```@eval
import MCMRSimulator.CLI: run_main_docs
run_main_docs("geometry create cylinders 1 geometry.json --radius 1 --repeats 2.5,2.5")
```

This will create a JSON file with the full information on the geometry:
```@eval
import Markdown
text = read("geometry.json", String)
Markdown.parse("```json\n$(text)\n```")
```

You can see how that the `repeats` and `radius` keywords have been set to our predefined values.
You can alter these and other geometry properties by editing this JSON directly or using the flags when creating the geometry.
For a full overview of these flags, you can run:
```
mcmr geometry create cylinders --help
```
```@eval
import MCMRSimulator.CLI: run_main_docs
run_main_docs("geometry create cylinders --help")
```
How these various properties affect the simulation is described [here](@ref properties).

The procedure to create [`Walls`](@ref), [`Spheres`](@ref), or [`Annuli`](@ref) is very similar as for the [`Cylinders`](@ref) illustrated above.
Randomly distributed cylinders, annuli, and spheres can be created using `mcmr geometry create-random`.

## Defining the sequence
In this case we will use an idealised DWI sequence produced by [`MRIBuilder.jl`](https://open.win.ox.ac.uk/pages/ndcn0236/mribuilder.jl/dev/).
Note that any [pulseq](http://pulseq.github.io/) file can be used instead.
```@example
using MRIBuilder
seq = DWI(TE=80, bval=2., Δ=40, δ=:min, TR=2000)
write_sequence("dwi.seq", seq)

# Optional to plot the sequence
using CairoMakie
f = plot_sequence(seq)
f
save("tutorial_cli_sequence.png", f); # hide
nothing # hide
```
![](tutorial_cli_sequence.png)

## Running the simulation
To get instructions on running the simulations, we can check the help message of `mcmr run`:
```bash
mcmr run --help
```
```@eval
import MCMRSimulator.CLI: run_main_docs
run_main_docs("run --help")
```

We can see that in addition to defining the geometry and the sequence, we can also control the simulation properties such as the `--diffusivity`, `--R1`, and `--R2`.

The simulation is initialised by randomly distributing a number of spins (set by `--Nspins`) uniformly across a bounding box with size given by `--voxel-size`.
This initial state might also contain bound spins (if the `--density` flag was set to a non-zero value during the geometry generation).

The DWI sequence defined above contains a `MRIBuilder.SingleReadout` object at the echo time (80 ms). By default, this is used for readout:
```bash
mcmr run geometry.json dwi.seq -o signal.csv
```
```@eval
import MCMRSimulator.CLI: run_main_docs
run_main_docs("run geometry.json dwi.seq -o signal.csv --seed=1")
```

This produces the CSV file, which looks like
```@eval
import Markdown
text = read("signal.csv", String)
Markdown.parse("```\n$(text)\n```")
```

The columns in this file store the following information:
- "sequence": integer; index of the sequence (always 1 if only single sequence used)
- "bvec": integer; index of the gradient orientation (if a `--bvec` flag is provided)
- "TR": integer; index of the repetition time that this data was acquired (between 1 and the value of `--nTR`)
- "readout": integer; index of the readout within a TR.
- "subset": integer; index of the subset of the total signal (e.g., intra-axonal) that has been output (see the `--subset` flag). The total signal will always be included with "subset" of 0.
- "nspins": total number of spins contributing to the signal (might change for certain subsets of spins)
- "longitudinal": average longitudinal signal
- "transverse": average transverse signal
- "phase": average phase of the signal (in degrees)
- "Sx": signal strength in the x-direction
- "Sy": signal strength in the y-direction

We can also output the signal of specific subsets of spins. For example, in the following we request to separately the output for just the spins inside the cylinders and just the spins outside of the cylinders.
```bash
mcmr run geometry.json dwi.seq -o signal.csv --subset inside --subset outside
```
```@eval
import MCMRSimulator.CLI: run_main_docs
run_main_docs("run geometry.json dwi.seq -o signal2.csv --subset inside --subset outside --seed=2")
```

We can see two additional rows in the output. 
These new rows are the in same order as the `--subset` flags provided to `mcmr run` and can be distinguished based on the "subset" column.
```@eval
import Markdown
text = read("signal2.csv", String)
Markdown.parse("```\n$(text)\n```")
```
We can see in the second row that inside the cylinders the transverse signal is very close to the number of spins, 
indicating that there has been very little dephasing due to the diffusion weighting inside the cylinders.
On the other hand, we did lose most of the signal outside of the cylinders (i.e., the transverse signal is much lower than the number of spins in the third row).
All the spins are either inside or outside the cylinders, so in this case the first row is simply the sum of the next two.

A more complete state of all the spins can be produced using the `--output-snapshot` flag.
For example, the command
```bash
mcmr run geometry.json dwi.seq --output-snapshot snapshot.csv
```
```@eval
import MCMRSimulator.CLI: run_main_docs
run_main_docs("run geometry.json dwi.seq --output-snapshot snapshot.csv --seed=3")
```
will produce a file named "snapshot.csv" with:
```@eval
import Markdown
lines = split(read("snapshot.csv", String), '\n')
text = join(lines[1:5], '\n')
Markdown.parse("```csv\n$(text)\n...\n```")
```
Each row corresponds to the state of a single spin. In addition to all the columns listed above, we now have 4 more columns:
- "spin": integer; index of the spin
- "x"/"y"/"z": floats; position of the spin at the time of the readout

The readout times can be adjusted using the `--nTR`, `--time`, and `--skip-TR` flags.
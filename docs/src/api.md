```@meta
CurrentModule = MCMRSimulator
```

# [MCMRSimulator.jl public API](@id api)
This is the public API for [MCMRSimulator](https://git.fmrib.ox.ac.uk/ndcn0236/MCMRSimulator.jl).
For a more user-friendly introduction, click [here](@ref Introduction).

## Running simulations
```@docs
Simulation
readout
evolve
Subset
Snapshot
BoundingBox
Spin
MCMRSimulator.get_subset
MCMRSimulator.SpinOrientation
MCMRSimulator.SpinOrientationSum
MCMRSimulator.TimeStep
```

## Creating geometry
### Geometry types
```@docs
Annuli
MCMRSimulator.Annulus
Cylinders
MCMRSimulator.Cylinder
Walls
MCMRSimulator.Wall
Spheres
MCMRSimulator.Sphere
Mesh
MCMRSimulator.Triangle
BendyCylinder
MCMRSimulator.Ring
```
### Geometry helper functions
```@docs
load_mesh
random_positions_radii
```

## Querying simulation output
```@docs
position
longitudinal
transverse
phase
orientation
get_sequence
isinside
stuck
stuck_to
get_time
```

## Plotting
```@docs
PlotPlane
plot_snapshot
plot_geometry
plot_trajectory
plot_off_resonance
simulator_movie
```

## Probing MCMR internals
```@docs
MCMRSimulator.gyromagnetic_ratio
MCMRSimulator.get_rotation
MCMRSimulator.get_readouts
MCMRSimulator.IndexedReadout
MCMRSimulator.susceptibility_off_resonance
MCMRSimulator.ObstructionGroup
MCMRSimulator.IndexedObstruction
MCMRSimulator.nvolumes
MCMRSimulator.fix
MCMRSimulator.fix_susceptibility
MCMRSimulator.FixedGeometry
MCMRSimulator.FixedObstructionGroup
MCMRSimulator.surface_relaxation
MCMRSimulator.surface_density
MCMRSimulator.dwell_time
MCMRSimulator.permeability
MCMRSimulator.GlobalProperties
MCMRSimulator.R1
MCMRSimulator.R2
MCMRSimulator.off_resonance
MCMRSimulator.install_cli
MCMRSimulator.run_main
```

## Internal modules
The internals of these sub-modules are considered private and might change at any time.

Each of these modules corresponds to a file in the gitlab repository

```@docs
MCMRSimulator
MCMRSimulator.Constants
MCMRSimulator.Methods
MCMRSimulator.Spins
MCMRSimulator.Evolve
MCMRSimulator.Relax
MCMRSimulator.Properties
MCMRSimulator.Plot
MCMRSimulator.Subsets
MCMRSimulator.SequenceParts
MCMRSimulator.Simulations
MCMRSimulator.TimeSteps
MCMRSimulator.CLI
MCMRSimulator.CLI.Geometry
MCMRSimulator.CLI.Run
MCMRSimulator.Geometries
MCMRSimulator.Geometries.User
MCMRSimulator.Geometries.User.Obstructions
MCMRSimulator.Geometries.User.Obstructions.Fields
MCMRSimulator.Geometries.User.Obstructions.ObstructionTypes
MCMRSimulator.Geometries.User.Obstructions.ObstructionGroups
MCMRSimulator.Geometries.User.RandomDistribution
MCMRSimulator.Geometries.User.Fix
MCMRSimulator.Geometries.User.SizeScales
MCMRSimulator.Geometries.User.FixSusceptibility
MCMRSimulator.Geometries.User.LoadMesh
MCMRSimulator.Geometries.User.JSON
MCMRSimulator.Geometries.User.ToMesh
MCMRSimulator.Geometries.User.SplitMesh
MCMRSimulator.Geometries.Internal
MCMRSimulator.Geometries.Internal.Properties
MCMRSimulator.Geometries.Internal.RayGridIntersection
MCMRSimulator.Geometries.Internal.BoundingBoxes
MCMRSimulator.Geometries.Internal.Intersections
MCMRSimulator.Geometries.Internal.Obstructions
MCMRSimulator.Geometries.Internal.Obstructions.FixedObstructions
MCMRSimulator.Geometries.Internal.Obstructions.ObstructionIntersections
MCMRSimulator.Geometries.Internal.Obstructions.Rounds
MCMRSimulator.Geometries.Internal.Obstructions.Shifts
MCMRSimulator.Geometries.Internal.Obstructions.Triangles
MCMRSimulator.Geometries.Internal.Obstructions.Walls
MCMRSimulator.Geometries.Internal.FixedObstructionGroups
MCMRSimulator.Geometries.Internal.Reflections
MCMRSimulator.Geometries.Internal.Susceptibility
MCMRSimulator.Geometries.Internal.Susceptibility.Base
MCMRSimulator.Geometries.Internal.Susceptibility.Cylinder
MCMRSimulator.Geometries.Internal.Susceptibility.Triangle
MCMRSimulator.Geometries.Internal.Susceptibility.Annulus
MCMRSimulator.Geometries.Internal.Susceptibility.Grid
MCMRSimulator.Geometries.Internal.HitGrids
```
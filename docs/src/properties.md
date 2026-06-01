# [Simulation properties](@id properties)
How the spins behave is determined by the tissue [geometry](@ref geometry), the applied MRI sequences, and user-provided flags determining how the spin magnetisation evolves. Here we discuss how the spin magnetisation evolution can be affected by these user-provided flags.

For example, one such flag is the `diffusivity`, which can be set as a keyword argument while generating the [`Simulation`](@ref).
## MRI properties
MRI properties determine the spin evolution for free and stuck particles. They include:
- the longitudinal relaxation rate `R1`
- the transverse relaxation rate `R2`
- the global `off_resonance` field (i.e., any off-resonance not caused by the sequence or the geometry)
At the [`Simulation`](@ref) level these parameters can be set by supplying the `R1`, `R2`, or `off_resonance` flags (see [`MCMRSimulator.GlobalProperties`](@ref)), such as:
```julia
simulation = Simulation(sequences, R2=1/80)
```
These MRI properties can be locally altered when defining the [geometry](@ref geometry). In the geometry they can be seperately set for spins stuck to the geometry surface or those spins that are inside specific objects in the geometry. The total relaxation rate (and off-resonance field) is set by the sum of the global value, the value set for any surface the spin is stuck to, and the value set for any obstruction that the spin is inside of. A single spin might be inside of multiple obstructions at once, if they overlap. In that case, all of the overlapping compartments will be considered. For the off-resonance field there might also be a contribution of the magnetic suscpetibility of any [`Cylinders`](@ref), [`Annuli`](@ref), or [`Mesh`](@ref).

If not set at the global or local level, there will be no longitudinal or transverse relaxation and there will be no off-resonance field.

From the command line interface, the global parameters are set during the `mcmr run` command using `--diffusivity`, `--R1`, and `--R2` keywords.
Local parameters will already have been set at an earlier stage during the creation of the geometry using `mcmr geometry create/create-random`.

## Collision properties
Collision properties determine the behaviour of spins at the time of a collision. Like MRI properties they can be set at the global level (while creating [`Simulation`](@ref)) or overwritten at the local level ([geometry](@ref geometry)). There are four such properties:
- [`MCMRSimulator.surface_relaxation`](@ref): the rate with which transverse signal is lost at every collision. This rate is multiplied by the square root of the timestep to ensure the actual attenuation is timestep-indeptendent.
- [`MCMRSimulator.permeability`](@ref): the rate of spins passing through the surface (arbitrary units). Set to infinity for a purely permeable surface and to zero for an impermeable surface (default). If the spins do not pass through, they will undergo regular reflection (or get stuck, see below). Like `MT_fraction` it will be adjusted to take into account the timestep.
- [`MCMRSimulator.surface_density`](@ref) and [`MCMRSimulator.dwell_time`](@ref): These control the density and dwell time of spins on the surface. Depending on the MRI properties assigned to these stuck particles (see above), these stuck particles can be used to represent water stuck at the membranes due to surface tension or spins in the membrane itself (which is in exchange with the free water through magnetisation transfer).


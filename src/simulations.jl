"""
Defines the main [`Simulation`](@ref) object.
"""
module Simulations
import StaticArrays: SVector, SizedVector
import MRIBuilder: Sequence
import ..Geometries: ObstructionGroup, fix, fix_susceptibility
import ..Geometries.Internal: FixedGeometry, FixedObstruction, FixedSusceptibility, susceptibility_off_resonance, prepare_isinside!
import ..Spins: Spin, Snapshot, SpinOrientation, stuck
import ..Methods: get_time
import ..Properties: GlobalProperties, R1, R2, off_resonance
import ..TimeSteps: TimeStep

"""
    Simulation(
        sequences; geometry=[], diffusivity=3.,
        R1=0, T1=Inf, R2=0, T2=Inf, off_resonance=0, MT_fraction=0, permeability=0,, 
        timestep=<default parameters>,
    )

Defines the setup of the simulation and stores the output of the run.

# Argument
## General parameters:
- `sequences`: Vector of `MRIBuilder.Sequence` objects. During the spin random walk the simulation will keep track of the spin magnetisations for all of the provided sequences.
- `geometry`: Set of obstructions, which can be used to restrict the diffusion, produce off-resonance fields, alter the local T1/T2 relaxation, and as sources of magnetisation transfer.
- `diffusivity`: Rate of the random motion of the spins in um^2/ms (default: 3).
- `verbose`: set to false to silence descriptions of the simulation parameters (default: true).

## MRI properties
These parameters determine the evolution and relaxation of the spin magnetisation.
- `R1`/`T1`: sets the longitudinal relaxation rate (R1 in kHz) or relaxation time (T1=1/R1 in ms). This determines how fast the longitudinal magnetisation returns to its equilibrium value of 1.
- `R2`/`T2`: sets the transverse relaxation rate (R2 in kHz) or relaxation time (T2=1/R2 in ms). This determines how fast the transverse magnetisation is lost.
- `off_resonance`: Size of the off-resonance field in this voxel in kHz.
These MRI properties can be overriden for spins inside the [`ObstructionGroup`](@ref) objects of the `geometry`.

## Collision parameters
These parameters determine how parameters behave when hitting the [`ObstructionGroup`](@ref) objects of the `geometry`.
They can be overriden for individual objects for each [`ObstructionGroup`](@ref).
- `MT_fraction`: the fraction of magnetisation transfered between the obstruction and the water spin at each collision.
- `permeability`: the rate of spins passing through the surface in arbitrary units (set to infinity for fully permeable surface).
- `surface_density`: Density of spins stuck on the surface relative to the volume density of hte free water.
- `dwell_time`: Typical time that spins will stay at the surface after getting stuck.
Note that `MT_fraction` and `permeability` are internally adjusted to make their effect independent of the timestep.

## Timestep parameters
`timestep` controls the timepoints at which the simulation is evaluated. 
By default, the maximum allowable timestep will be determined by the geometry and biophysical parameters as described in [`MCMRSimulator.TimeStep`](@ref MCMRSimulator.TimeSteps.TimeStep).
That documentation also describes how to adjust these settings.
The timestep can also be set to a number to ignore any of these parameters.
Note that a too large timestep will lead to inaccurate results.

# Running the simulation
To run a [`Snapshot`](@ref) of spins through the simulations you can use one of the following functions:
- `evolve`: evolves the spins in the snapshot until a single given time and returns that state in a new [`Snapshot`](@ref).
- `readout`: evolves the spins to particular times in each TR and return the total signal at that time (or a [`Snapshot`](@ref)).
"""
struct Simulation{N, NG, G<:FixedGeometry{NG}, IG<:FixedGeometry, O<:FixedSusceptibility}
    # N sequences, datatype T
    sequences :: SizedVector{N, Sequence}
    diffusivity :: Float64
    properties :: GlobalProperties
    geometry :: G
    inside_geometry :: IG
    susceptibility :: O
    timestep :: TimeStep
    flatten::Bool
    verbose::Bool
    function Simulation(
        sequences, 
        diffusivity::Float64,
        properties::GlobalProperties,
        geometry::FixedGeometry,
        inside_geometry::FixedGeometry,
        susceptibility::FixedSusceptibility,
        timestep::TimeStep,
        flatten::Bool,
        verbose::Bool
    )
        nseq = length(sequences)

        new{nseq, length(geometry), typeof(geometry), typeof(inside_geometry), typeof(susceptibility)}(
            SizedVector{nseq, Sequence}(sequences),
            diffusivity,
            properties,
            geometry,
            inside_geometry,
            susceptibility,
            timestep,
            flatten,
            verbose,
        )
    end
end

function Simulation(
    sequences;
    diffusivity=3.,
    geometry=[],
    permeability=0.,
    surface_density=0.,
    dwell_time=0.,
    surface_relaxation=0.,
    R1=0.,
    R2=0.,
    off_resonance=0.,
    verbose=true,
    timestep=(),
)
    flatten = false
    if isa(sequences, Sequence)
        sequences = [sequences]
        flatten = true
    elseif length(sequences) == 0
        sequences = Sequence[]
    end
    susceptibility = fix_susceptibility(geometry)
    geometry = fix(geometry; permeability=permeability, density=surface_density, dwell_time=dwell_time, relaxation=surface_relaxation)

    inside_geometry = filter(geometry) do obstruction
        ~all(all(iszero.(getproperty(obstruction.volume, s))) for s in (:R1, :R2, :off_resonance))
    end
    prepare_isinside!.(inside_geometry)

    default_properties = GlobalProperties(; R1=R1, R2=R2, off_resonance=off_resonance)
    if iszero(diffusivity) && length(geometry) > 0
        @warn "Restrictive geometry will have no effect, because the diffusivity is set at zero"
    end
    return Simulation(
        sequences, 
        Float64(diffusivity),
        default_properties,
        geometry,
        inside_geometry,
        susceptibility,
        timestep isa Number ? TimeStep(timestep, Inf) : TimeStep(; verbose=verbose, diffusivity=diffusivity, geometry=geometry, timestep...),
        flatten,
        verbose,
    )
end

function Base.show(io::IO, sim::Simulation{N}) where {N}
    function print_geometry()
        print(io, "Geometry(")
        if length(sim.geometry) > 20
            print(io, "$(length(sim.geometry)) groups of obstructions")
        else
            for obstruction in sim.geometry
                print(io, string(obstruction), ", ")
            end
        end
        print(io, ")")
    end
    if get(io, :compact, false)
        seq_text = sim.flatten ? "single sequence" : "$N sequences"
        print(io, "Simulation($seq_text, ")
        print_geometry()
        print(io, ", D=$(sim.diffusivity)um^2/ms, $(sim.properties))")
    else
        print(io, "Simulation(")
        print_geometry()
        print(io, ", D=$(sim.diffusivity)um^2/ms, $(sim.properties)):\n")
        print(io, "$N sequences:\n")
        for seq in sim.sequences
            print(io, seq)
        end
    end
end

function Snapshot(nspins::Integer, simulation::Simulation{N}, bounding_box=500; kwargs...) where {N}
    Snapshot(nspins, bounding_box, simulation.geometry; nsequences=N, kwargs...)
end
_to_snapshot(spins::Int, simulation::Simulation, bounding_box) = Snapshot(spins, simulation, bounding_box)
_to_snapshot(spins::AbstractVector{<:Real}, simulation::Simulation, bounding_box) = _to_snapshot(Spin(position=spins), simulation, bounding_box)
_to_snapshot(spins::AbstractVector{<:AbstractVector{<:Real}}, simulation::Simulation, bounding_box) = _to_snapshot([Spin(position=pos) for pos in spins], simulation, bounding_box)
function _to_snapshot(spins::AbstractMatrix{<:Real}, simulation::Simulation, bounding_box) 
    if size(spins, 2) != 3
        spins = transpose(spins)
    end
    @assert size(spins, 2) == 3
    _to_snapshot([spins[i, :] for i in 1:size(spins, 1)], simulation, bounding_box)
end
_to_snapshot(spins::Spin, simulation::Simulation, bounding_box) = _to_snapshot([spins], simulation, bounding_box)
_to_snapshot(spins::AbstractVector{<:Spin}, simulation::Simulation, bounding_box) = _to_snapshot(Snapshot(spins), simulation, bounding_box)
_to_snapshot(spins::Snapshot{1}, simulation::Simulation{nseq}, bounding_box) where {nseq} = nseq == 1 ? deepcopy(spins) : Snapshot(spins, nseq)
_to_snapshot(spins::Snapshot{N}, simulation::Simulation{N}, bounding_box) where {N} = deepcopy(spins)

produces_off_resonance(sim::Simulation) = produces_off_resonance(sim.geometry)
propose_times(sim::Simulation, t_start, t_end) = propose_times(sim.time_controller, t_start, t_end, sim.sequences, sim.diffusivity)

for symbol in (:R1, :R2, :off_resonance)
    @eval begin
        function $symbol(spin_or_pos, simulation::Simulation)
            $symbol(spin_or_pos, simulation.geometry, simulation.properties)
        end
    end
end

"""
    susceptibility_off_resonance(simulation, spin)

Computes the susceptibility off-resonance caused by all susceptibility sources in the [`Simulation`](@ref) affecting the [`Spin`](@ref)

The field is computed in ppm. Knowledge of the scanner `B0` is needed to convert it into KHz.
"""
susceptibility_off_resonance(simulation::Simulation, spin::Spin) = susceptibility_off_resonance(simulation, spin.position, stuck(spin) ? spin.reflection.inside : nothing)
susceptibility_off_resonance(simulation::Simulation, position::AbstractVector, inside::Union{Nothing, Bool}=nothing) = susceptibility_off_resonance(simulation.susceptibility, SVector{3, Float64}(position), inside)
function susceptibility_off_resonance(simulation::Simulation, old_pos::AbstractVector, new_pos::AbstractVector) 
    if iszero(length(simulation.susceptibility))
        return 0.
    end
    along_tract = rand()
    susceptibility_off_resonance(simulation.susceptibility, along_tract * old_pos .+ (1 - along_tract) .* new_pos, nothing)
end


end
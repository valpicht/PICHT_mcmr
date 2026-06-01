"""
Defines the functions that run the actual simulation:
- [`readout`](@ref): get total signal or [`Snapshot`](@ref) at any `MRIBuilder.ADC` objects in the sequences.
- [`evolve`](@ref): Return a single [`Snapshot`](@ref) with the state of the simulation at a given time. This snapshot can be used as initialisation for further runs.

All of these functions call `evolve_to_time` under the hood to actually run the simulation.
"""
module Evolve
import StaticArrays: SVector, MVector
import LinearAlgebra: norm, ⋅
import MRIBuilder: Sequence, variables, B0, build_sequence
import MRIBuilder.Components: InstantGradient, InstantPulse
import Rotations
import Bessels: besseli0
import ..SequenceParts: SequencePart, MultSequencePart, InstantSequencePart, iter_parts, get_readouts, IndexedReadout, first_TR_with_all_readouts, NoGradient
import ..Methods: get_time
import ..Spins: @spin_rng, Spin, Snapshot, stuck, SpinOrientationSum, get_sequence, orientation, SpinOrientation, static_vector_type
import ..Simulations: Simulation, _to_snapshot
import ..Relax: relax!
import ..Properties: GlobalProperties, stick_probability
import ..Subsets: Subset, get_subset
import ..Geometries.Internal: Reflection, detect_intersection, empty_intersection, has_intersection, surface_relaxation, permeability, surface_density, direction, previous_hit, dwell_time, empty_reflection, FixedGeometry

"""
Supertype for any Readout accumulator.
"""
abstract type ReadoutAccumulator end
abstract type SingleAccumulator <: ReadoutAccumulator end

"""
    TotalSignalAccumulator(readout_time)

Accumulator used to store the total signal.

The actual state of the `Snapshot` is discarded (unlike [`SnapshotAccumulator`](@ref)),
which can be used by setting `return_snapshot=true` in [`readout`](@ref).

After the simulation a [`SpinOrientationSum`](@ref) will be returned.
"""
struct TotalSignalAccumulator <: SingleAccumulator
    nspins :: Ref{Int}
    as_vector :: MVector{3, Float64}
    nwrite :: Ref{Int}
    TotalSignalAccumulator(time::Number) = new(
        Ref(0),
        zero(MVector{3, Float64}),
        Ref(0)
    )
end

"""
    SnapshotAccumulator(readout_time)

Accumulator used to store the full spin state.

This uses a lot of memory.
If you are just interested in the total signal use [`TotalSignalAccumulator`](@ref) instead
by setting `return_snapshot=false` in [`readout`](@ref).

After the simulation a [`Snapshot`](@ref) will be returned.
"""
struct SnapshotAccumulator{N, ST} <: SingleAccumulator
    spins :: Vector{Vector{Spin{N, ST}}}
    time :: Float64
    nwrite :: Ref{Int}
    function SnapshotAccumulator{N}(time::Number) where {N}
        st = static_vector_type(N){SpinOrientation}
        return new{N, st}(
            Vector{Spin{N, st}}[],
            Float64(time),
            Ref(0)
        )
    end
end

"""
    FillerAccumulator()

An accumulator that will be used to fill any gaps in the regular grid of the readout accumulation.

After the simulation `nothing` is returned.
"""
struct FillerAccumulator <: SingleAccumulator
end

"""
    GridAccumulator(simulation, start_time; kwargs...)

Creates a grid of accumulators to readout the results of the `simulation` starting at `start_time` ms.

Keywords match the relevant ones described in [`readout`](@ref).

# Properties
- `grid`: 4-dimensional grid of [`SingleAccumulator`](@ref) objects. Each dimension corresponds to a possible dimension of the [`readout`](@ref) output, namely:
    1. multiple sequences
    2. multiple readouts in single TR
    3. subsequent sequence repetitions
    4. individual subsets (see [`Subset`](@ref))
- `flatten`: length-4 vector of booleans indicating which dimensions should be kept when producing the output from [`readout`](@ref).
"""
struct GridAccumulator{T<:SingleAccumulator, S<:Simulation}
    grid :: Array{T, 4}
    subsets :: Vector{Subset}
    simulation :: S
    flatten :: SVector{4, Bool}
    last_readout :: Tuple{Int, Int, Int}
    first_TR :: Int
end

function GridAccumulator(simulation::Simulation{N}, start_time::Number; noflatten=false, subset=Subset(), return_snapshot=false, readouts=nothing, nTR=nothing, kwargs...) where {N}
    flatten_subset = subset isa Subset
    if flatten_subset
        subset = [subset]
    end

    if iszero(N)
        if isnothing(readouts)
            error("Should set `readouts` explicitly when running a simulation with 0 sequences.")
        end
        if !return_snapshot
            error("`return_snapshot` should be set to true when running a simulation with 0 sequences")
        end
        flatten_readouts = readouts isa Number
        use_readouts = flatten_readouts ? [readouts] : readouts

        grid = Array{SnapshotAccumulator{0}}(undef, 1, length(readouts), 1, length(subset))
        for i_r in 1:length(use_readouts)
            for i_s in 1:length(subset)
                grid[1, i_r, 1, i_s] = SnapshotAccumulator{0}(use_readouts[i_r])
            end
        end
        to_flatten = SVector{4, Bool}(
            noflatten ? 
            (false, false, false, false) :
            (
                true,
                flatten_readouts,
                true,
                flatten_subset,
            )
        )
        return GridAccumulator(
            grid,
            subset,
            simulation,
            to_flatten,
            (1, argmax(use_readouts), 1),
            0
        )
    end

    actual_readouts = collect.(get_readouts.(simulation.sequences, start_time; readouts=readouts, nTR=nTR, kwargs...))

    if any(iszero.(length.(actual_readouts)))
        error("No readouts scheduled for at least one of the sequences.")
    end
    flat_readouts = vcat(actual_readouts...)

    nreadouts = maximum(flat_readouts) do ro
        ro.readout
    end
    flatten_readouts = (
        isnothing(readouts) ?
        nreadouts == 1 :
        readouts isa Number
    )

    first_TR = minimum(flat_readouts) do ro
        ro.TR
    end - 1

    grid_size = (
        N,
        nreadouts,
        maximum(flat_readouts) do ro
            ro.TR
        end - first_TR,
        length(subset)
    )

    as_dict = Dict(
        (sequence, ro.readout, ro.TR) => ro
        for sequence in 1:N for ro in actual_readouts[sequence]
    )
    acc_type = return_snapshot ? SnapshotAccumulator{1} : TotalSignalAccumulator

    grid = Array{SingleAccumulator}(undef, grid_size...)
    for index in eachindex(IndexCartesian(), grid)
        i_seq, i_readout, i_TR, _ = Tuple(index)
        if (i_seq, i_readout, i_TR + first_TR) in keys(as_dict)
            grid[index] = acc_type(as_dict[(i_seq, i_readout, i_TR + first_TR)].time)
        else
            grid[index] = FillerAccumulator()
        end
    end
    if !any(isa.(grid, FillerAccumulator))
        grid = Array{acc_type}(grid)
    end
    to_flatten = SVector{4, Bool}(
        noflatten ? 
        (false, false, false, false) :
        (
            simulation.flatten,
            flatten_readouts,
            isnothing(nTR),
            flatten_subset,
        )
    )


    max_readout_time = maximum(flat_readouts) do ro
        ro.time
    end
    last_readout = [key for (key, value) in pairs(as_dict) if value.time == max_readout_time][1]
    return GridAccumulator(
        grid,
        subset,
        simulation,
        to_flatten,
        last_readout,
        first_TR,
    )
end

"""
    readout!(single_accumulator, spins)

Adds the `spins` to the readout of `single_accumulator`
"""
function readout!(acc::TotalSignalAccumulator, spins::Vector{<:Spin{1}})
    acc.nspins[] += length(spins)
    acc.as_vector .+= orientation(Snapshot(spins))
    acc.nwrite[] += 1
end

function readout!(acc::SnapshotAccumulator{N}, spins::Vector{<:Spin{N}}) where {N}
    push!(acc.spins, deepcopy(spins))
    acc.nwrite[] += 1
end

readout!(::FillerAccumulator, spins::Vector{<:Spin{1}}) = nothing

"""
    readout!(grid_accumulator, spins, index_sequence, readout)

Adds the `spins` to the accumulators corresponding to each subset for given sequence and readout.

Returns true if this is the last readout that needs to be considered in the simulation.
"""
function readout!(acc::GridAccumulator, spins::Vector{<:Spin{N}}, index_sequence::Int, readout::IndexedReadout) where {N}
    @assert N <= 1
    if readout.TR <= acc.first_TR
        return false
    end
    for (index_subset, subset) in enumerate(acc.subsets)
        use_spins = get_subset(spins, acc.simulation, subset)
        full_index = CartesianIndex(
            index_sequence,
            readout.readout,
            readout.TR - acc.first_TR,
            index_subset
        )
        readout!(acc.grid[full_index], use_spins)
    end
    return (index_sequence, readout.readout, readout.TR) == acc.last_readout
end


"""
    has_nwrite(accumulator, target)

Returns true if the number of times that each of the accumulators has been written matches `target`.
"""
has_nwrite(acc::SingleAccumulator, nwrite::Int) = acc.nwrite[] == nwrite
has_nwrite(acc::GridAccumulator, nwrite::Int) = all(has_nwrite.(acc.grid, nwrite))
has_nwrite(::FillerAccumulator, nwrite::Int) = true

"""
    fix_accumulator(single_accumulator)

Return the accumulated readout results as a `SpinOrientationSum` or `Snapshot`.
"""
fix_accumulator(acc::TotalSignalAccumulator) = SpinOrientationSum(
    SpinOrientation(acc.as_vector),
    acc.nspins[]
)
fix_accumulator(acc::SnapshotAccumulator) = Snapshot(vcat(acc.spins...), acc.time)
fix_accumulator(::FillerAccumulator) = nothing

"""
    fix_accumulator(grid_accumulator)

Return the accumulated readout results as an array.
"""
function fix_accumulator(acc::GridAccumulator)
    full_grid = fix_accumulator.(acc.grid)
    indices = map(acc.flatten) do flatten
        flatten ? 1 : (:)
    end
    return full_grid[indices...]
end


"""
    readout(spins, simulation[, readout_times]; bounding_box=<1x1x1 mm box>, skip_TR=0, nTR=1, return_snapshot=false, subset=<all>)

Evolves a set of spins through the [`Simulation`](@ref).
Returns the total signal or a full [`Snapshot`](@ref) at every readout time in the simulated sequences over one or more repetition times (TRs).

# Positional arguments:
- `spins`: Number of spins to simulate or an already existing [`Snapshot`](@ref).
- `simulation`: [`Simulation`](@ref) object defining the environment, scanner, and sequence(s).
- `times` (optional): time of the readouts relative to the start of the TR (in ms). If not provided, the times of any `MRIBuilder.ADC` objects in the sequence will be used (see [`get_readouts`](@ref) for details).

# Keyword arguments:
- `bounding_box`: size of the voxel in which the spins are initiated in um (default is 1000, corresponding to a 1x1x1 mm box centered on zero). Can be set to a [`BoundingBox`](@ref MCMRSimulator.Geometries.Internal.BoundingBoxes.BoundingBox) object for more control.
- `skip_TR`: Number of repetition times to skip before starting the readout. 
    Even if set to zero (the default), the simulator will still skip the current TR before starting the readout 
    if the starting snapshot is from a time past one of the sequence readouts.
    See [`get_readouts`](@ref) for details.
- `nTR`: number of TRs for which to store the output. See [`get_readouts`](@ref) for details.
- `return_snapshot`: set to true to output the state of all the spins as a [`Snapshot`](@ref) at each readout instead of a [`SpinOrientationSum`](@ref) with the total signal.
- `subset`: Return the signal/snapshot for a subset of all spins. Can be set to a single or a vector of [`Subset`](@ref) objects. If set to a vector, this will add an attional dimension to the output.

# Returns
The function returns an up to 3-dimensional (KxLxMxN) array, with the following dimensions:
- `K`: the number of sequences. This dimension is not included if the simulation only contains a single sequencen (and this single sequence is not passed into the [`Simulation`](@ref) as a vector).
- `L`: the number of readout times with a single TR. This dimension is skipped if the `readout_times` is set to a scalar number. This dimension might contain `nothing`s for sequences that contain fewer `Readout.ADC` objects than the maximum (`M`).
- `M`: the number of TRs (controlled by the `nTR` keyword). If `nTR` is not explicitly set by the user, this dimension is skipped.
- `N`: the number of subsets (controlled by the `subset` keyword). If `subset` is set to a single value (<all> by default), this dimension is skipped.
By default each element of this matrix is either a [`SpinOrientationSum`](@ref) with the total signal.
If `return_snapshot=true` is set, each element is the full [`Snapshot`](@ref) instead.
"""
readout(spins, simulation::Simulation, new_readout_times=nothing; bounding_box=500, kwargs...) = readout_internal(_to_snapshot(spins, simulation, bounding_box), simulation, new_readout_times; kwargs...)

function readout_internal(snapshot::Snapshot{N}, simulation::Simulation{N}, new_readout_times=nothing; kwargs...) where {N}
    if :readouts in keys(kwargs)
        error("readout timings should be set as the 3rd positional argument, not a keyword argument.")
    end
    accumulator = GridAccumulator(simulation, snapshot.time; readouts=new_readout_times, kwargs...)
    repeat = :nTR in keys(kwargs) || :skip_TR in keys(kwargs)
    run_readout!(snapshot, simulation, accumulator, Val(repeat); readouts=new_readout_times)
    return fix_accumulator(accumulator)
end


# Special case when `spins` is an integer value
# Only run a limited number of spins at a time to save memory
function readout(spins::Integer, simulation::Simulation{N}, new_readout_times=nothing; bounding_box=500, kwargs...) where {N}
    if :readouts in keys(kwargs)
        error("readout timings should be set as the 3rd positional argument, not a keyword argument.")
    end
    if iszero(N)
        nruns = 1
    else
        total_magnetisations = spins * N
        nruns = Int(div(total_magnetisations, 1e7, RoundUp))
        if nruns > spins
            nruns = spins
        end
        if simulation.verbose && nruns > 1
            @info "The $spins simulated spins will be split over $nruns independent runs to save memory."
        end
    end

    nspins_min = Int(div(spins, nruns, RoundDown))
    nruns_extra = spins - nspins_min * nruns

    accumulator = GridAccumulator(simulation, 0.; readouts=new_readout_times, kwargs...)
    repeat = :nTR in keys(kwargs) || :skip_TR in keys(kwargs)
    for nspins_run in [
        fill(nspins_min + 1, nruns_extra)...,
        fill(nspins_min, nruns - nruns_extra)...
    ]
        run_readout!(_to_snapshot(nspins_run, simulation, bounding_box), simulation, accumulator, Val(repeat); readouts=new_readout_times)
    end
    return fix_accumulator(accumulator)
end


"""
    run_readout!(snapshot, simulation, accumulators, repeat; kwargs...)

Runs the simulation starting from `snapshot` and filling the `accumulators`.
"""
function run_readout!(snapshot::Snapshot{N}, simulation::Simulation{N}, accumulator::GridAccumulator, repeat::Val; kwargs...) where {N}
    @assert N > 0
    B0s = map(B0, simulation.sequences)
    for part in iter_parts(simulation.sequences, snapshot.time, repeat, simulation.timestep; kwargs...)
        if process_sequence_step!(snapshot.spins, simulation, part, B0s, accumulator)
            return
        end
    end
    error("Simulation ended before all readout accumulators were filled.")
end

function run_readout!(snapshot::Snapshot{0}, simulation::Simulation{0}, accumulator::GridAccumulator, repeat::Val{R}; readouts, kwargs...) where {R}
    if R
        error("Cannot set `nTR` or `skip_TR` when there are no sequences.")
    end
    empty_sequence = build_sequence() do 
        Sequence([1.]) 
    end
    B0s = zero(SVector{0, Float64})
    for part in iter_parts([empty_sequence], snapshot.time, Val(false), simulation.timestep; readouts=readouts)
        if part isa MultSequencePart
            draw_step!(snapshot.spins, simulation, MultSequencePart(part.duration, NoGradient[]), B0s)
        else
            if readout!(accumulator, snapshot.spins, 1, part.instants[1])
                return
            end
        end
    end
    error("Simulation ended before all readout accumulators were filled.")
end



"""
    evolve(snapshot, simulation[, new_time]; bounding_box=<1x1x1 mm box>)

Evolves the [`Snapshot`](@ref) through the [`Simulation`](@ref) to a new time.
Returns a [`Snapshot`](@ref) at the new time, which can be used as a basis for further simulation.
"""
function evolve(spins, simulation::Simulation{N}, new_time; TR=nothing, bounding_box=500) where {N}
    snapshot = _to_snapshot(spins, simulation, bounding_box)
    if isnothing(TR)
        if snapshot.time > new_time
            error("Cannot evolve snapshot back in time.")
        end
        res = readout_internal(snapshot, simulation, new_time, return_snapshot=true)
    else
        if iszero(N)
            error("Cannot set `TR` in evolve without having a sequence in simulation.")
        end
        rep_time = variables.TR(simulation.sequences[1])
        if !all(variables.duration(s) ≈ rep_time for s in simulation.sequences)
            error("Cannot evolve snapshot for multiple repetition times, because the simulation contains sequences with different TRs. Please set a `new_time` explicitly.")
        end
        if new_time > rep_time
            if isapprox(new_time, rep_time, rtol=1e-3)
                new_time = rep_time
            else
                error("When evolving to a time within a given `TR`, the time cannot be greater than the sequence repetition time.")
            end
        end
        current_TR = first_TR_with_all_readouts(simulation.sequences[1], snapshot.time; readouts=new_time)
        res = readout_internal(snapshot, simulation, new_time; skip_TR=TR - current_TR, return_snapshot=true)
    end
    if simulation.flatten || iszero(N)
        return res
    end
    @assert size(res) == (N, )
    spins = [
        Spin(
            res[1].spins[index].position,
            SVector{N, SpinOrientation}([
                r.spins[index].orientations[1]
                for r in res
            ]),
            res[1].spins[index].reflection,
            res[1].spins[index].rng,
        )
        for index in 1:length(snapshot.spins)
    ]
    return Snapshot(
        spins,
        res[1].time
    )
end

function process_sequence_step!(spins, simulation, sequence_part::MultSequencePart, B0s, _) 
    draw_step!(spins, simulation, sequence_part, B0s)
    return false
end
process_sequence_step!(spins, simulation, sequence_part::InstantSequencePart, _, accumulator) = apply_instants!(spins, sequence_part, accumulator)

"""
    draw_step!(spin(s), simulation, mult_sequence_part, B0s)

Updates the spin based on a random movement through the given geometry for a given `timestep`:
- draws the next location of the particle after `timestep` with given `simulation.diffusivity`.  
  This displacement will take into account the obstructions in `simulation.geometry`.
- The spin orientation will be affected by relaxation (see [`relax!`](@ref)) and potentially by magnetisation transfer during collisions.
"""
function draw_step!(spins::Vector{<:Spin{N}}, simulation::Simulation{N}, sequence_part::MultSequencePart{N}, B0s::AbstractVector{Float64}) where {N}
    Threads.@threads for spin in spins
        draw_step!(spin, simulation, sequence_part, B0s)
    end
end


function draw_step!(spin::Spin{N}, simulation::Simulation{N}, parts::MultSequencePart{N}, B0s::AbstractVector{Float64}, test_new_pos=nothing) where {N}
    if ~isnothing(test_new_pos)
        all_positions = [spin.position]
    end
    timestep = parts.duration
    if iszero(timestep)
        if ~isnothing(test_new_pos)
            return all_positions
        else
            return
        end
    end
    current_pos = spin.position
    is_stuck = stuck(spin)
    fraction_timestep = 0.

    found_solution = false
    @spin_rng spin begin
        if !stuck(spin)
            if isnothing(test_new_pos)
                rand_base_vec = randn(SVector{3, Float64})
                displacement = rand_base_vec .* sqrt(2 * simulation.diffusivity * timestep)
                new_pos = current_pos + displacement
            else
                new_pos = SVector{3, Float64}(test_new_pos)
            end
            reflection = Reflection(norm(new_pos - current_pos) / sqrt(2 * simulation.diffusivity * timestep))
            phit = (0, 0, false)
        end
        for _ in 1:1000000
            if is_stuck
                td = dwell_time(simulation.geometry, spin.reflection)
                fraction_stuck = -log(rand()) * td / timestep
                relax!(spin, spin.reflection.inside, simulation, parts, fraction_timestep, min(one(Float64), fraction_timestep + fraction_stuck), B0s)
                fraction_timestep += fraction_stuck
                if fraction_timestep >= 1
                    found_solution = true
                    break
                end
                phit = previous_hit(spin.reflection)
                reflection = spin.reflection
                displacement = direction(reflection, (1 - fraction_timestep) * timestep, simulation.diffusivity)
                new_pos = current_pos .+ displacement

                spin.reflection = empty_reflection
                is_stuck = false
            end

            collision = detect_intersection(
                simulation.geometry,
                current_pos,
                new_pos,
                phit,
            )

            use_distance = collision === empty_intersection ? 1. : max(prevfloat(collision.distance), 0.)
            if collision === empty_intersection
                next_fraction_timestep = 1.
                collision_pos = new_pos
            else
                next_fraction_timestep = fraction_timestep + (1 - fraction_timestep) * use_distance
                collision_pos = map((p1, p2) -> use_distance * p1 + (1 - use_distance) * p2, new_pos, current_pos)
            end

            # spin relaxation
            relax!(spin, collision_pos, simulation, parts, fraction_timestep, next_fraction_timestep, B0s)

            if ~has_intersection(collision)
                spin.position = new_pos
                found_solution = true
                break
            end

            relaxation = surface_relaxation(simulation.geometry, collision)
            if ~iszero(relaxation)
                collision_attenuation = exp(-sqrt(timestep) * relaxation)
                for orientation in spin.orientations
                    orientation.transverse *= collision_attenuation
                end
            end

            normed_rate = sqrt(timestep) * permeability(simulation.geometry, collision)
            permeability_prob = isinf(normed_rate) ? 1. : (1. - exp(-normed_rate) * besseli0(normed_rate))
            passes_through = isone(permeability_prob) || !(iszero(permeability_prob) || rand() > permeability_prob)
            reflection = Reflection(collision, new_pos - current_pos, reflection.ratio_displaced, 
                reflection.time_moved + (1 - fraction_timestep) * use_distance * timestep, 
                reflection.distance_moved + norm(new_pos - current_pos) * use_distance, 
                passes_through
            )
            current_pos = spin.position = collision_pos
            if ~isnothing(test_new_pos)
                push!(all_positions, current_pos)
            end

            sd = surface_density(simulation.geometry, collision)
            if !iszero(sd) && rand() < stick_probability(sd, dwell_time(simulation.geometry, collision), simulation.diffusivity, timestep)
                spin.reflection = reflection
                is_stuck = true
            else
                new_pos = current_pos .+ direction(reflection, (1 - next_fraction_timestep) * timestep, simulation.diffusivity)
            end

            phit = previous_hit(reflection)
            fraction_timestep = next_fraction_timestep
        end
    end
    if !found_solution
        error("Bounced single particle for 1000000 times in single step; terminating!")
    end
    if ~isnothing(test_new_pos)
        push!(all_positions, spin.position)
        return all_positions
    end
end

"""
    apply_instants!(spins, instants, accumulator)

Apply a set of `N` instants to a vector of spins for each of the `N` sequences being simulated.

Each instant can be:
- `nothing`: do nothing
- `MRIBuilder.Components.InstantPulse`: apply RF pulse rotation
- `MRIBuilder.Components.InstantGradient`: add phase corresponding to gradient
- `IndexedReadout`: add the current snapshot to `accumulator`

Returns `true` if this is the final readout and we should stop.
"""
apply_instants!(spins::Vector{<:Spin{N}}, instants::InstantSequencePart{N}, accumulator::GridAccumulator) where {N} = apply_instants!(spins, instants.instants, accumulator)
function apply_instants!(spins::Vector{<:Spin{N}}, instants::AbstractVector, accumulator::GridAccumulator) where {N}
    final = false
    for i in 1:N
        final |= apply_instants!(spins, i, instants[i], accumulator)
    end
    return final
end

apply_instants!(spins::Vector{<:Spin{N}}, instants::AbstractVector{Nothing}, _) where {N} = false
apply_instants!(spins::Vector{<:Spin}, index::Int, ::Nothing, _) = false

function apply_instants!(spins::Vector{<:Spin}, index::Int, grad::InstantGradient, _)
    Threads.@threads for spin in spins
        new_phase = rad2deg(spin.position ⋅ variables.qvec(grad))
        spin.orientations[index].phase += new_phase
    end
    return false
end

function apply_instants!(spins::Vector{<:Spin}, index::Int, readout::IndexedReadout, accumulator::GridAccumulator)
    readout!(accumulator, get_sequence.(spins, index), index, readout)
end

function apply_instants!(spins::Vector{<:Spin{1}}, index::Int, readout::IndexedReadout, accumulator::GridAccumulator)
    readout!(accumulator, spins, index, readout)
end

function apply_instants!(spins::Vector{<:Spin}, index::Int, pulse::InstantPulse, _)
    rotation = Rotations.RotationVec(
        deg2rad(pulse.flip_angle) * cosd(pulse.phase),
        deg2rad(pulse.flip_angle) * sind(pulse.phase),
        0.
    )
    Threads.@threads for spin in spins
        orient = spin.orientations[index]
        new_orient = SpinOrientation(rotation * orientation(orient))
        orient.longitudinal = new_orient.longitudinal
        orient.transverse = new_orient.transverse
        orient.phase = new_orient.phase
    end
    return false
end

end
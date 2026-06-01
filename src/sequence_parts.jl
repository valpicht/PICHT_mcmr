"""
Breaks up one or more `MRIBuilder.Sequence` objects into individual timesteps.

The main interface is by calling `iter_parts`, which will result into a vector of some subtypes of `SequencePart`.

To just get the readouts call [`MCMRSimulator.get_readouts`](@ref MCMRSimulator.SequenceParts.get_readouts).
"""
module SequenceParts
import StaticArrays: SVector
import LinearAlgebra: norm
import MRIBuilder: BaseSequence, BaseBuildingBlock, waveform_sequence, events, get_gradient, edge_times, get_pulse, iter_instant_gradients, iter_instant_pulses, make_generic, variables, Wait
import MRIBuilder.Components: GradientWaveform, RFPulseComponent, NoGradient, ConstantGradient, ChangingGradient, InstantGradient, InstantPulse, split_timestep, EventComponent, SingleReadout
import ..TimeSteps: TimeStep
import ..Spins: static_vector_type


"""
A short part of the sequence that can be handled by the simulator.

There are 4 types:
- [`EmptyPart`](@ref): no pulse or gradient
- [`ContantPart`](@ref): no pulse; constant gradient
- [`LinearPart`](@ref): no pulse; linear gradient
- [`PulsePart`](@ref): RF pulse; constant gradient

To break down a generic sequence in these parts will require some approximations.
"""
abstract type SequencePart end

abstract type NoPulsePart <: SequencePart end

struct EmptyPart <: NoPulsePart
end


struct ConstantPart <: NoPulsePart
    strength :: SVector{3, Float64}
end


struct LinearPart <: NoPulsePart
    start :: SVector{3, Float64}
    final :: SVector{3, Float64}
end


"""
    ConstantPulse(amplitude, phase, frequency)

Stores the RF pulse state during a single timestep.

It contains:
- `amplitude` of RF pulse in kHz (assumed to be constant over the timestep).
- `phase` of the RF pulse at the beginning of the timestep in degrees.
- off-resonance `frequency` of the RF pulse in kHz.
"""
struct ConstantPulse
    amplitude :: Float64
    phase :: Float64
    frequency :: Float64
end

struct PulsePart{T<:NoPulsePart} <: SequencePart
    pulse :: Vector{ConstantPulse}
    gradient :: T
end


"""
    MultSequencePart(duration, parts)

A set of N [`SequencePart`](@ref) objects representing overlapping parts of `N` sequences.
"""
struct MultSequencePart{N, T<:SequencePart, ST<:AbstractVector{T}}
    duration :: Float64
    parts :: ST
    function MultSequencePart(duration::Number, parts::AbstractVector)
        N = length(parts)
        if iszero(N)
            return new{0, EmptyPart, SVector{0, EmptyPart}}(Float64(duration), SVector{0, EmptyPart}())
        end

        base_type = typeof(parts[1])
        if all(p -> p isa base_type, parts)
            T = base_type
        else
            T = SequencePart
        end
        return new{N, T, static_vector_type(N){T}}(
            Float64(duration),
            static_vector_type(N){T}(parts)
        )
    end
end

"""
    InstantSequencePart(instants)

A set of `N` instant pulses/gradients that should be applied to the spins.

Some of the instants might be `nothing`.
"""
struct InstantSequencePart{N, T, ST<:AbstractVector{T}}
    instants :: ST
    function InstantSequencePart(instants::AbstractVector)
        T = Union{typeof.(instants)...}
        N = length(instants)
        return new{N, T, static_vector_type(N){T}}(static_vector_type(N){T}(instants))
    end
end

"""
    IndexedReadout(time, TR, readout)

Represents a readout in an MR sequence.

Generate these for a sequence of interest using [`get_readouts`](@ref).

All indices are integers. They refer to:
- `time`: time since beginning of sequence in ms
- `TR`: which TR the simulation is in (defaults to 0 if not set).
- `readout`: which readout within the TR (or total sequence) this is.
"""
struct IndexedReadout
    time::Float64
    TR::Int
    readout::Int
end


iter_building_blocks_raw(seq::BaseSequence) = Iterators.flatten(iter_building_blocks_raw.(seq))
iter_building_blocks_raw(bb::BaseBuildingBlock) = [bb]

iter_building_blocks_no_time(seq, repeat::Val{false}) = Iterators.flatten([iter_building_blocks_raw(seq), [Wait(Inf)]])
iter_building_blocks_no_time(seq, repeat::Val{true}) = Iterators.cycle(Iterators.flatten([iter_building_blocks_raw(seq), [:TR]]))

iter_building_blocks(seq::BaseSequence, repeat) = Iterators.accumulate(iter_building_blocks_no_time(seq, repeat); init=(1, 0., nothing, Ref{Int}(0))) do (prev_TR, prev_time, prev_bb, readout_index), bb
    if isnothing(prev_bb)
        return (prev_TR, prev_time, bb, Ref{Int}(0))
    elseif bb == :TR
        return (prev_TR + 1, 0., Wait(0.), Ref{Int}(0))
    else
        return (prev_TR, prev_time + variables.duration(prev_bb), bb, readout_index)
    end
end


function iter_part_times(seq::BaseSequence, repeat::Val; readouts=nothing)
    rep_time = variables.TR(seq)
    Iterators.flatten(
        Iterators.map(iter_building_blocks(seq, repeat)) do (TR, time, bb, readout_index)
            gradients = Tuple{Float64, Float64, GradientWaveform}[]
            pulses = Tuple{Float64, Float64, RFPulseComponent}[]
            instants = Tuple{Float64, Union{InstantGradient, InstantPulse, SingleReadout, IndexedReadout}}[]

            full_time = (TR - 1) * rep_time + time

            current_grad = NoGradient(0.)
            current_time = time
            next_time = time
            for key in keys(bb)
                component = bb[key]
                if component isa GradientWaveform
                    duration = component.duration
                    current_time = next_time
                    next_time += duration
                    if duration > 0.
                        push!(gradients, (current_time, next_time, component))
                    end
                elseif component isa Tuple{<:Number, <:EventComponent}
                    delay, event = component
                    if event isa Union{InstantGradient, InstantPulse}
                        push!(instants, (current_time + delay, event))
                    elseif event isa RFPulseComponent
                        duration = variables.duration(event)
                        push!(pulses, (current_time + delay, current_time + delay + duration, event))
                    elseif event isa SingleReadout && isnothing(readouts)
                        readout_index[] += 1
                        push!(instants, (current_time + delay, IndexedReadout((TR - 1) * rep_time + current_time + delay, TR, readout_index[])))
                    end
                end
            end

            if !isnothing(readouts)
                # add custom readouts
                end_time = time + variables.duration(bb)
                for (index, readout) in enumerate(readouts)
                    if (
                        (time < readout <= end_time) ||  # readout is during this block
                        (iszero(time) && iszero(readout) && time != end_time) # or readout is at beginning of TR (and this block is not instant)
                    )
                        push!(instants, (readout, IndexedReadout((TR - 1) * rep_time + readout, TR, index)))
                    end
                end
            end

            et = sort!(unique!([
                time,
                [t for (_, t, _) in gradients]...,
                [t for (t, _, _) in pulses]...,
                [t for (_, t, _) in pulses]...,
                [t for (t, _) in instants]...,
            ]))

            current_grad(time) = findfirst(gradients) do (t1, t2, _)
                t2 > time
            end
            current_pulse(time) = findfirst(pulses) do (t1, t2, _)
                t1 <= time < t2
            end

            result = Tuple{Int, Float64, Float64, Tuple{Float64, GradientWaveform}, Union{Nothing, Tuple{Float64, RFPulseComponent}}, Union{Nothing, InstantGradient, InstantPulse, SingleReadout, IndexedReadout}}[]
            for (t1, t2) in zip(et[1:end-1], et[2:end])
                index_grad = current_grad(t1)
                t_grad, _, grad = gradients[index_grad]
                index_pulse = current_pulse(t1)
                pulse = isnothing(index_pulse) ? nothing : (t1 - pulses[index_pulse][1], pulses[index_pulse][end])
                for (t, instant) in instants
                    if t == t1
                        push!(result, (TR, t1, t1, (t1 - t_grad, grad), pulse, instant))
                    end
                end
                push!(result, (TR, t1, t2, (t1 - t_grad, grad), pulse, nothing))
            end
            for (t, instant) in instants
                if t == et[end]
                    t_grad, _, grad = iszero(length(gradients)) ? (t, t, NoGradient(0.)) : gradients[end]
                    push!(result, (TR, t, t, (t - t_grad, grad), nothing, instant))
                end
            end
            return result
        end
    )
end

struct _IterEdges{R, T}
    individual_iterators::Vector{T}
    TRs::Vector{Float64}
    tstart::Float64
end

Base.iterate(ie::_IterEdges) = Base.iterate(ie, (ie.tstart, iterate.(ie.individual_iterators)))
Base.IteratorSize(::Type{<:_IterEdges{false}}) = Base.SizeUnknown()
Base.IteratorSize(::Type{<:_IterEdges{true}}) = Base.IsInfinite()

function Base.iterate(ie::_IterEdges, my_state)
    current_time, states = my_state

    new_states = map(ie.individual_iterators, ie.TRs, states) do iter, rep_time, state
        (iTR, _, t2, _, _, instant), _ = state
        while isnothing(instant) && ((t2 + (iTR - 1) * rep_time) <= current_time || isapprox(t2, current_time, atol=1e-8))
            state = Base.iterate(iter, state[2])
            (iTR, _, t2, _, _, instant), _ = state
        end
        return state
    end
    if any(!isnothing(state[1][end]) for state in new_states)
        result = (map(new_states) do state
            instant = state[1][end]
            return instant
        end)
        new_states = map(ie.individual_iterators, new_states) do iter, state
            if isnothing(state[1][end])
                return state
            else
                return Base.iterate(iter, state[2])
            end
        end
        return (result, (current_time, new_states))
    end

    next_time = minimum(zip(new_states, ie.TRs); init=Inf) do (((iTR, _, t2, _, _, _), _), rep_time)
        t2 + (iTR - 1) * rep_time
    end
    if isinf(next_time)
        return nothing
    end
    return ((current_time, next_time, map(new_states, ie.TRs) do state, rep_time
        (iTR, t1, _, (t_grad, grad), pulse, _), _ = state
        correct_pulse = isnothing(pulse) ? pulse : (current_time - t1 - (iTR - 1) * rep_time + pulse[1], pulse[2])
        return ((current_time - t1 - (iTR - 1) * rep_time + t_grad, grad), correct_pulse)
    end), (next_time, new_states))
end

function iter_part_times(sequences::AbstractVector{<:BaseSequence}, tstart::Number, repeat::Val; kwargs...)
    iter_part_times(collect(sequences), tstart, repeat; kwargs...)
end

function iter_part_times(sequences::Vector{<:BaseSequence}, tstart::Number, repeat::Val{R}; kwargs...) where {R}
    iters = iter_part_times.(sequences, repeat; kwargs...)
    TRs = Float64.(variables.TR.(sequences))
    dropped = [Iterators.dropwhile(i) do (iTR, _, t2, _, _, _)
        ((iTR - 1) * TR + t2) < tstart
    end for (i, TR) in zip(iters, TRs)]
    return _IterEdges{R, eltype(dropped)}(
        dropped,
        TRs,
        Float64(tstart)
    )
end

function iter_part_times(sequences::AbstractVector{<:BaseSequence}, tstart, repeat::Val, timestep::TimeStep; kwargs...)
    Iterators.flatten(
        Iterators.map(iter_part_times(sequences, tstart, repeat; kwargs...)) do var
            if !(var isa Tuple)
                # instants
                return (var, )
            end
            
            (t1, t2, bbs) = var

            max_grad = iszero(length(sequences)) ? 0. : maximum(bbs) do ((t_grad, grad), _)
                max(
                    norm(variables.gradient_strength(grad, t_grad)),
                    norm(variables.gradient_strength(grad, t_grad + t2 - t1))
                )
            end
            use_timestep = timestep(max_grad)
            nsteps = isinf(use_timestep) ? 1 : Int(div(t2 - t1, use_timestep, RoundUp))
            splits = range(t1, t2, length=nsteps+1)
            return (
                (
                    t1p, t2p, map(bbs) do ((t_grad, grad), pulse)
                        correct_pulse = isnothing(pulse) ? nothing : (pulse[1] + t1p - t1, pulse[2])
                        return ((t_grad + t1p - t1, grad), correct_pulse)
                    end
                )
                for (t1p, t2p) in zip(splits[1:end-1], splits[2:end])
            )
        end
    )
end


function iter_parts(sequences::AbstractVector{<:BaseSequence}, tstart::Number, repeat::Val, timestep::TimeStep; kwargs...)
    Iterators.map(iter_part_times(sequences, tstart, repeat, timestep; kwargs...)) do var
        if !(var isa Tuple)
            return InstantSequencePart(var)
        else
            (t1, t2, bbs) = var
            if iszero(length(sequences))
                return MultSequencePart(t2 - t1, SVector{0, EmptyPart}())
            end
            return MultSequencePart(t2 - t1, map(bbs) do ((t_grad, gradient), gp)
                if gradient isa NoGradient
                    grad_part = EmptyPart()
                elseif gradient isa ConstantGradient
                    grad_part = ConstantPart(variables.gradient_strength(gradient))
                elseif gradient isa ChangingGradient
                    grad_part = LinearPart(variables.gradient_strength(gradient, t_grad), variables.gradient_strength(gradient, t_grad + t2 - t1))
                else
                    error("Gradient waveform $gradient is not implemented in the MCMR simulator yet.")
                end
                if isnothing(gp)
                    return grad_part
                else
                    (pulse_t1, pulse) = gp
                    pulse_t2 = (t2 - t1) + pulse_t1
        
                    max_ts = split_timestep(pulse, 1e-3)
        
                    nparts = isinf(max_ts) ? 1 : Int(div(pulse_t2 - pulse_t1, max_ts, RoundUp))
                    ts = (pulse_t2 - pulse_t1) / nparts
        
                    time_parts = ((1:nparts) .- 0.5) .* ts .+ pulse_t1
                    parts = [
                        ConstantPulse(variables.amplitude(pulse, t), variables.phase(pulse, t) - variables.frequency(pulse, t) * 180 * ts, variables.frequency(pulse, t))
                        for t in time_parts
                    ]
                    return PulsePart{typeof(grad_part)}(parts, grad_part)
                end
            end)
        end
    end
end

nreadouts_per_TR(seq::BaseSequence) = length(collect(Iterators.filter(iter_part_times(seq, Val(false))) do state
    state[end] isa SingleReadout
end))

function first_TR_with_all_readouts(seq::BaseSequence, start_time::Number; readouts=nothing)
    if iszero(start_time)
        return 1
    else
        current_TR = 1
        for ro in get_readouts(seq, 0.; nTR=typemax(Int), readouts=readouts)
            if ro.time <= start_time
                current_TR = ro.TR + 1
            else
                return current_TR
            end
        end
    end
    error()
end

"""
    get_readouts(sequence, start_time; readouts=nothing, nTR=1, skip_TR=0)

Returns a iterator of the readouts ([`IndexedReadout`](@ref) objects) that will be used for the given sequence in the simulator.

This can be used to identify which readouts will be (or have been) used in the simulation by running:
`collect(get_readouts(sequence, snapshot.current_time; kwargs...))`
where `snapshot` is the starting snapshot (which has a `current_time` of 0 by default) and `kwargs` are the keyword arguments used in
[`readout`](@ref MCMRSimulator.Evolve.readout) (i.e., `readouts, `nTR`, and `skip_TR`).

By default the readout/ADC objects with the actual sequence definition are used.
These can be overriden by `readouts`, which can be set to a vector of the timings of the readouts within each TR.

# Non-repeating sequences
This is the behaviour if both `nTR` and `skip_TR` are not set by the user.
Any readouts before the `start_time` are ignored.
Any readouts after the `start_time` (whether from the `readouts` keyword or within the sequence definition) are returned.

# Repeating sequences
This is the behaviour if either `nTR` or `skip_TR` or both are not set by the user.
If at `start_time` any of the readouts in the current TR have already passed, then the readout will only start in the next TR (unless `skip_TR` is set to -1).
We will skip an additional number of TRs given by `skip_TR` (default: 0).
Then readouts will continue for the number of TRs given by `nTR` (default: 1).
"""
function get_readouts(seq::BaseSequence, start_time::Number; readouts=nothing, nTR=nothing, skip_TR=nothing) 
    repeat = !(isnothing(nTR) && isnothing(skip_TR))
    if isnothing(nTR)
        nTR = 1
    end
    if isnothing(skip_TR)
        skip_TR = 0
    end

    if !isnothing(readouts)
        readouts = collect(readouts)
        if repeat
            rep_time = variables.TR(seq)
            max_ro = maximum(readouts)
            if max_ro > rep_time
                if isapprox(max_ro, rep_time, rtol=1e-3)
                    @warn "Adjusting maximum readout time ($max_ro) to match the TR of $rep_time."
                else
                    error("Readouts have been scheduled at $max_ro ms beyond the sequence repetitition time of $rep_time ms.")
                end
            end
            readouts = min.(readouts, rep_time)
        end
    end
    base_iterator = iter_part_times(seq, Val(repeat); readouts=readouts)
    if repeat
        first_TR = first_TR_with_all_readouts(seq, start_time, readouts=readouts) + skip_TR
        last_TR = first_TR + nTR - 1
        drop_first = Iterators.dropwhile(res -> res[1] < first_TR, base_iterator)
        drop_last = Iterators.takewhile(res -> res[1] <= last_TR, drop_first)
    else
        drop_last = Iterators.dropwhile(res -> res[2] < start_time, base_iterator)
    end

    filtered = Iterators.filter(drop_last) do res
        return res[end] isa IndexedReadout
    end
    return Iterators.map(filtered) do res
        return res[end]
    end
end


end
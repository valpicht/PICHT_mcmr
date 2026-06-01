"""
Implements the Bloch equations to model the spin magnetisation under influence of a sequence during a single timestep.
"""
module Relax

import StaticArrays: SVector
import Rotations
import LinearAlgebra: ⋅
import ..Constants: gyromagnetic_ratio
import ..Spins: Spin, SpinOrientation, stuck, R1, R2, off_resonance, stuck_to, orientation
import ..Geometries.Internal: susceptibility_off_resonance, MRIProperties
import ..Properties: GlobalProperties
import ..SequenceParts: MultSequencePart, SequencePart, NoPulsePart, EmptyPart, ConstantPart, LinearPart, PulsePart, ConstantPulse
import ..Simulations: Simulation


const NewPosType = Union{SVector{3, Float64}, Nothing, Bool}

"""
    relax!(spin, new_pos, simulations, sequence_parts, t1, t2, B0s)

Evolves a `spin` during part of a sequence represented a [`MultSequencePart`](@ref) (containing one [`SequencePart`](@ref) for each [`SpinOrientation`](@ref) in the [`Spin`](@ref)).

The spin is evolved for the duration from `t1` and `t2`, where `t=0` corresponds to the start of the sequence part and `t=1` to the end.
The spin is assumed to travel from `spin.position` to `new_pos` during this period.

The spin will precess around an effective RF pulse and relax with a given `R1` and `R2`.
The R1 and R2 are determined based on the [`Simulation`](@ref) parameters, potentially influenced by the geometry.
The effective RF pulse is determined by the sequence and any changes in the off-resonance field due to the geometry or set in the [`Simulation`](@ref).
"""
function relax!(spin::Spin{N}, new_pos::NewPosType, simulation::Simulation{N}, parts::MultSequencePart{N}, t1::Float64, t2::Float64, B0s::AbstractVector{Float64}) where {N}

    mean_pos = (spin.position .+ new_pos) ./ 2
    props = MRIProperties(simulation.geometry, simulation.inside_geometry, simulation.properties, mean_pos, spin.reflection)

    off_resonance_unscaled = susceptibility_off_resonance(simulation, spin.position, new_pos) * 1e-6 * gyromagnetic_ratio

    # pre-compute the R1 and R2 attenuation.
    # This will be applied for any sequences for which there is no RF pulse active.
    # Pre-computing this can greatly speed up runs with large number of sequences.
    pre_comp = (
        R1_att=exp(-props.R1 * (t2 - t1) * parts.duration),
        R2_att=exp(-props.R2 * (t2 - t1) * parts.duration),
    )

    for index in 1:N
        relax!(spin.orientations[index], spin.position, new_pos, parts.parts[index], props, parts.duration, t1, t2, off_resonance_unscaled * B0s[index], pre_comp)
    end
end


# no active RF pulse
function relax!(orient::SpinOrientation, old_pos::SVector{3, Float64}, new_pos::NewPosType, part::NoPulsePart, props::MRIProperties, duration::Float64, t1::Float64, t2::Float64, off_resonance::Float64, pre_comp::NamedTuple)
    full_off_resonance = off_resonance + props.off_resonance + grad_off_resonance(part, old_pos, new_pos, t1, t2)
    timestep = (t2 - t1) * duration
    orient.phase += full_off_resonance * timestep * 360
    relax_single_step!(orient, props, pre_comp)
end

function relax_single_step!(orient::SpinOrientation, props::MRIProperties, pre_comp::NamedTuple)
    orient.longitudinal = (1 - (1 - orient.longitudinal) * pre_comp.R1_att) 
    orient.transverse *= pre_comp.R2_att
end

function relax_single_step!(orient::SpinOrientation, props::MRIProperties, timestep::Float64)
    orient.longitudinal = (1 - (1 - orient.longitudinal) * exp(-props.R1 * timestep)) 
    orient.transverse *= exp(-props.R2 * timestep)
end

# with active RF pulse
function relax!(orient::SpinOrientation, old_pos::SVector{3, Float64}, new_pos::NewPosType, pulse::PulsePart, props::MRIProperties, duration::Float64, t1::Float64, t2::Float64, off_resonance::Float64, ::NamedTuple)
    relax_time = 1 / max(props.R1, props.R2)
    if isinf(relax_time)
        nsplit_rotation = Val(1)
    else
        internal_timestep = 1/length(pulse.pulse)
        nsplit_rotation = Val(Int(div(internal_timestep * duration, relax_time/10, RoundUp)))
    end
    relax!(orient, old_pos, new_pos, pulse, props, duration, t1, t2, off_resonance, nsplit_rotation)
end
    
function relax!(orient::SpinOrientation, old_pos::SVector{3, Float64}, new_pos::NewPosType, pulse::PulsePart, props::MRIProperties, duration::Float64, t1::Float64, t2::Float64, off_resonance::Float64, split_rotation::Val)
    started = iszero(t1)
    internal_timestep = 1/length(pulse.pulse)
    single_pulse_duration = duration * internal_timestep
    call_apply_pulse!(part, index, t1_use, t2_use) = apply_pulse!(
        orient, part, props, single_pulse_duration, t1_use / internal_timestep - index + 1, t2_use / internal_timestep - index + 1,
        off_resonance + props.off_resonance + grad_off_resonance(pulse.gradient, old_pos, new_pos, t1_use, t2_use), split_rotation
    )

    for (index, part) in enumerate(pulse.pulse)
        t_int = internal_timestep * index
        if !started
            if t_int > t1
                if t_int > t2
                    call_apply_pulse!(part, index, t1, t2)
                    return
                else
                    call_apply_pulse!(part, index, t1, t_int)
                end
                started = true
            end
        else
            if t_int > t2
                call_apply_pulse!(part, index, t_int - internal_timestep, t2)
                return
            else
                call_apply_pulse!(part, index, t_int - internal_timestep, t_int)
            end
        end
    end
    @assert isone(t2)
end

function apply_pulse!(orient::SpinOrientation, pulse::ConstantPulse, props::MRIProperties, duration::Float64, t1::Float64, t2::Float64, full_off_resonance::Float64, ::Val{1})
    flip_angle = pulse.amplitude * 2π * duration * (t2 - t1)

    if iszero(t1)
        start_phase = pulse.phase
    else
        start_phase = pulse.phase + pulse.frequency * duration * 360 * t1
    end

    rotation = Rotations.RotationVec(
        flip_angle * cosd(start_phase),
        flip_angle * sind(start_phase),
        (full_off_resonance - pulse.frequency) * 2π * duration * (t2 - t1)
    )

    relax_single_step!(orient, props, duration * (t2 - t1) / 2)
    new_orient = SpinOrientation(rotation * orientation(orient))
    orient.longitudinal = new_orient.longitudinal
    orient.transverse = new_orient.transverse
    orient.phase = new_orient.phase + (pulse.frequency * 360 * duration * (t2 - t1))
    relax_single_step!(orient, props, duration * (t2 - t1) / 2)
end

function apply_pulse!(orient::SpinOrientation, pulse::ConstantPulse, props::MRIProperties, duration::Float64, t1::Float64, t2::Float64, full_off_resonance::Float64, ::Val{N}) where {N}
    # relaxation times are short compared with rotation
    internal_stepsize = (t2 - t1) / N
    for i in 1:N
        t2_sub = t1 + i * internal_stepsize
        t1_sub = t2_sub - internal_stepsize
        apply_pulse!(orient, pulse, props, duration, t1_sub, t2_sub, full_off_resonance, Val(1))
    end
end

grad_off_resonance(grad::EmptyPart, old_pos::SVector{3, Float64}, new_pos::NewPosType, old_time::Float64, new_time::Float64) = 0.

grad_off_resonance(grad::ConstantPart, old_pos::SVector{3, Float64}, new_pos::SVector{3, Float64}, old_time::Float64, new_time::Float64) = ((old_pos .+ new_pos) ⋅ grad.strength) / 2
grad_off_resonance(grad::ConstantPart, old_pos::SVector{3, Float64}, ::Union{Nothing, Bool}, old_time::Float64, new_time::Float64) = old_pos ⋅ grad.strength

function grad_off_resonance(grad::LinearPart, old_pos::SVector{3, Float64}, new_pos::SVector{3, Float64}, old_time::Float64, new_time::Float64)
    old_grad = iszero(old_time) ? grad.start : (old_time .* grad.final .+ (1 - old_time) .* grad.start)
    new_grad = isone(new_time) ? grad.final : (new_time .* grad.final .+ (1 - new_time) .* grad.start)
    return (
        (2 .* old_pos .+ new_pos) ⋅ old_grad +
        (old_pos .+ 2 .* new_pos) ⋅ new_grad
    ) / 6
end

function grad_off_resonance(grad::LinearPart, old_pos::SVector{3, Float64}, ::Union{Nothing, Bool}, old_time::Float64, new_time::Float64)
    mean_time = (old_time + new_time) / 2
    mean_grad = (mean_time .* grad.final .+ (1 - mean_time) .* grad.start)
    return mean_grad ⋅ old_pos
end

end
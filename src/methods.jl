"""
Defines methods shared across multiple sub-modules.

- [`get_time`](@ref)
- `norm_angle`
- [`MCMRSimulator.get_rotation`](@ref MCMRSimulator.Methods.get_rotation)
- `off_resonance`
"""
module Methods

import Rotations
import LinearAlgebra: norm, cross, ⋅, I
import StaticArrays: SMatrix

"""
    get_time(snapshot)
    get_time(sequence_component)
    get_time(sequence, sequence_index)

Returns the time in milliseconds that a snapshot was taken or that a sequence component will have effect.
"""
function get_time end


"""
    norm_angle(angle)

Normalises an angle in degrees, so that it is between it is in the range (-180, 180]
"""
function norm_angle(angle)
    angle = mod(angle, 360)
    if angle > 180
        angle -= 360
    end
    angle
end


function off_resonance end

"""
    get_rotation(rotation_mat, ndim)

Returns a (3, `ndim`) rotation matrix, that is the relevant part of the full 3x3 `rotation_mat` to map to the x-axis (if `ndim` is 1) or the x-y plance (if `ndim` is 2).
If `ndim` is 3, the full rotation matrix `rotation_mat` is returned.

    get_rotation(vector, ndim; reference_dimension)

Returns the (3, `ndim`) rotation matrix mapping the `vector` to the `reference_dimension`.
By default, the `reference_dimension is the x-direction (if `ndim` is 1 or 3) or the z-direction (if `ndim` is 2).
`vector` and `reference_dimension` can be a length 3 array or one of the symbols :x, :y, or :z (representing vectors in those cardinal directions).
"""
get_rotation(rotation::Rotations.Rotation, ndim::Int; reference_dimension=nothing) = get_rotation(Rotations.RotMatrix(rotation), ndim)
get_rotation(rotation::Rotations.RotMatrix, ndim::Int; reference_dimension=nothing) = get_rotation(rotation.mat, ndim)
function get_rotation(rotation::AbstractMatrix, ndim::Int; reference_dimension=nothing)
    if size(rotation) == (3, 3) && ndim < 3
        rotation = rotation[:, 1:ndim]
    end
    SMatrix{3, ndim, Float64}(rotation)
end

function get_rotation(rotation::AbstractVector, ndim::Int; reference_dimension=nothing)
    T = typeof(rotation[1])
    return get_rotation(T.(rotation), ndim)
end

function get_rotation(rotation::AbstractVector{<:AbstractVector}, ndim::Int; reference_dimension=nothing)
    return get_rotation(hcat(rotation...), ndim)
end

function get_rotation(rotation::AbstractVector{<:Number}, ndim::Int; reference_dimension=nothing)
    @assert length(rotation)==3
    nr = norm(rotation)
    if iszero(nr)
        error("Cannot rotate to vector of [0, 0, 0]")
    end
    normed = rotation / nr
    if ndim == 1
        return get_rotation(reshape(normed, 3, 1), 1)
    end
    try_vec = [0., 1., 0.]
    vec1 = cross(normed, try_vec)
    if iszero(norm(vec1))
        try_vec = [1., 0, 0]
        vec1 = cross(normed, try_vec)
    end
    if isnothing(reference_dimension)
        reference = dimension_symbol_to_vec(:I, ndim)
    elseif isa(reference_dimension, Symbol)
        reference = dimension_symbol_to_vec(reference_dimension, ndim)
    else
        reference = reference_dimension
    end

    # rotate reference to normed
    angle = acos(normed ⋅ reference)
    vec_rotation = cross(reference, normed)
    if iszero(norm(vec_rotation))
        return get_rotation(abs(angle) < 1 ? I(3) : -I(3), ndim)
    end
    return get_rotation(Rotations.AngleAxis(angle, vec_rotation...), ndim)
end

function get_rotation(rotation::Symbol, ndim::Int; reference_dimension=nothing)
    return get_rotation(dimension_symbol_to_vec(rotation, ndim), ndim; reference_dimension=reference_dimension)
end

function dimension_symbol_to_vec(s::Symbol, ndim)
    return Dict(
        :x => [1, 0, 0],
        :y => [0, 1, 0],
        :z => [0, 0, 1],
        :I => (ndim == 2 ? [0, 0, 1] : [1, 0, 0]),
    )[s]
end

end
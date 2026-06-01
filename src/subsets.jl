"""
Support for selecting a subset of the total signal.

Types:
- [`Subset`](@ref)

Functions:
- [`get_subset`](@ref)
"""
module Subsets
import Base: @kwdef
import ..Spins: Snapshot, Spin, stuck_to
import ..Simulations: Simulation
import ..Geometries.Internal: FixedGeometry, isinside

_arguments = """
- `bound`: set to true to return only bound spins, to false to return only free spins (default: whether spins are bound is not relevant).
- `inside`: set to true to return only spins inside the geometry, to false to return only spins outside of the geometry (default: whether spins are inside or outside is not relevant). It can also be set to a positive integer number. Only spins that are inside that exact number of obstructions will be returned.
- `geometry_index`: set to an index to only consider that obstruction group within the total geometry (default: consider the full geometry).
- `obstruction_index`: set to an index to only consider that obstruction group within the total geometry (default: consider the full geometry).
"""

"""
    Subset(; bound=nothing, inside=nothing, geometry_index=nothing, obstruction_index=nothing)

This creates a helper object to extract a subset of a [`Snapshot`](@ref) from the total snapshot.
It defines which spins should be kept.
This definition is determined by:
$(_arguments)
"""
@kwdef struct Subset
    bound :: Union{Nothing, Bool} = nothing
    inside :: Union{Nothing, Bool, Int} = nothing
    geometry_index :: Union{Nothing, Int} = nothing
    obstruction_index :: Union{Nothing, Int} = nothing
end
    

"""
    get_subset(snapshot, simulation, subset)
    get_subset(snapshot, simulation; bound=nothing, inside=nothing, geometry_index=nothing, obstruction_index=nothing)

Returns a subset of the [`Snapshot`](@ref) from the [`Simulation`](@ref) that obey some specific properties.
These properties can be either defined by a [`Subset`](@ref) object or a set of keyword arguments.

These keyword arguments are:
$(_arguments)
"""
get_subset(spins::AbstractVector{<:Spin}, args...; kwargs...) = get_subset(Snapshot(spins), args...; kwargs...).spins

get_subset(snapshot::Snapshot, simulation::Union{Simulation, FixedGeometry}; kwargs...) = get_subset(snapshot, simulation, Subset(; kwargs...))

get_subset(snapshot::Snapshot, simulation::Simulation, subset::Subset) = get_subset(snapshot, simulation.geometry, subset)

function get_subset(snapshot::Snapshot{N}, geometry::FixedGeometry, subset::Subset) where {N}
    if isnothing(subset.geometry_index) && !isnothing(subset.obstruction_index)
        error("Cannot select a specific obstruction from an obstruction group, if no obstruction group has been selected. Please, do not set `obsruction_index`, without also setting `geometry_index`.")
    end
    if isnothing(subset.bound) && isnothing(subset.inside)
        return snapshot
    end
    if ~isnothing(subset.geometry_index)
        geometry = filter(og -> og.original_index == subset.geometry_index, geometry)
    end
    if iszero(length(geometry))
        if subset.inside in (nothing, false, 0) && subset.bound in (nothing, false)
            return snapshot
        else
            return Snapshot(0, time=snapshot.time, nsequences=N)
        end
    end
    include = ones(Bool, length(snapshot))
    fixed_geometry_indices = [g.parent_index for g in geometry]
    if ~isnothing(subset.bound)
        function _isbound(spin::Spin)
            (g_index, o_index) = stuck_to(spin)
            return g_index in fixed_geometry_indices && (isnothing(subset.obstruction_index) || o_index == subset.obstruction_index)
        end
        include .&= (_isbound.(snapshot) .== subset.bound)
    end
    if ~isnothing(subset.inside)
        function _number_isinside(spin::Spin)
            nmatch = 0
            for g in geometry
                if !isnothing(subset.geometry_index) && subset.geometry_index != g.original_index
                    continue
                end
                o_indices = isinside(g, spin.position, spin.reflection)
                valid_indices = filter(o_indices) do index
                    isnothing(subset.obstruction_index) || subset.obstruction_index == index
                end
                nmatch += length(valid_indices)
            end
            return nmatch
        end
        _isinside(spin::Spin) = _number_isinside(spin) > 0
        if subset.inside isa Bool
            include .&= (_isinside.(snapshot) .== subset.inside)
        else
            include .&= (_number_isinside.(snapshot) .== subset.inside)
        end
    end
    return snapshot[include]
end


end

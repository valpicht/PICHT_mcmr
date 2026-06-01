"""
Defines parent type and functions for all magnetic susceptibility sources.


Types:
- `BaseSusceptibility`

Functions:
- `single_susceptibility`
- `single_susceptibility_gradient`
"""
module Base

"""
    BaseSusceptibility{N}

Parent type of all susceptibility sources.
Each `BaseSusceptibility` object represents a single susceptibility source in `N`-dimensional space 
(e.g., N=2 for [`CylinderSusceptibility`](@ref) or N=3 for [`PointSusceptibility`](@ref)).
"""
abstract type BaseSusceptibility{N} end

"""
    single_susceptibility(source::BaseSusceptibility, offset, distance, stuck_to, b0_field)

Computes the off-resonance field contribution from a [`BaseSusceptibility`](@ref) source given:
- `offset` vector from spin position to the centre of the `source`
- `distance` of the spin from the source (which is the norm of `offset`)
- `stuck_to` `nothing` for free particle. Otherwise, boolean indicating whether the particle is stuck to the inside or outside of the obstructions it is stuck to.
- `b0_field` N-length vector with the magnetic field within the coordinate system of the susceptibility source
"""
function single_susceptibility end

"""
    single_susceptibility_gradient(source::BaseSusceptibility)

Maximal off-resonance gradient in ppm/um produced by a single susceptibility source.
"""
function single_susceptibility_gradient end


end

"""
Off-resonance field due to magnetic susceptibility of an hollow cylinder, which is approximated as being infinitely thin.
"""
module Cylinder
import StaticArrays: SVector
import ..Base: BaseSusceptibility, single_susceptibility, single_susceptibility_gradient

"""
    CylinderSusceptibility(radius, g_ratio, chi_I, chi_A, b0_field)

Creates a cylindrical susceptibility source, which produces the intra- and extra-axonal field expected for a myelin sheath.
More realistic myelin sheaths with a finite width can be modeled using `AnnulusSusceptibility`.
"""
struct CylinderSusceptibility <: BaseSusceptibility{2}
    rsq :: Float64
    internal_field :: Float64
    external_field :: Float64
end


function CylinderSusceptibility(radius::Number, g_ratio::Number, chi_I::Number, chi_A::Number, b0_field::AbstractVector)
    @assert length(b0_field) == 2
    sin_theta_sq = sum(b0_field .* b0_field)
    CylinderSusceptibility(
        radius^2,
        -0.75 * chi_A * log(g_ratio) * sin_theta_sq,
        2 * (chi_I + chi_A/4) * (1 - g_ratio^2) / (1 + g_ratio)^2 * radius^2 * sin_theta_sq,
    )
end

"""
    single_susceptibility(cylinder, position, distance, stuck_to, b0_field)

Computed by the hollow cylinder fiber model from [Wharton_2012](@cite).
The myelin sheath is presumed to be an infinitely thin cylinder.
"""
function single_susceptibility(cylinder::CylinderSusceptibility, position::AbstractVector, distance::Number, stuck_inside::Union{Nothing, Bool}, b0_field::SVector{2, Float64})
    if iszero(cylinder.internal_field) && iszero(cylinder.external_field)
        return zero(Float64)
    end
    rsq = distance * distance
    if ~isnothing(stuck_inside) && rsq ≈ cylinder.rsq
        inside = stuck_inside
    else
        inside = rsq < cylinder.rsq
    end
    if inside
        return cylinder.internal_field
    else
        sin_theta_sq = b0_field[1] * b0_field[1] + b0_field[2] * b0_field[2]
        cos2 = (b0_field[1] * position[1] + b0_field[2] * position[2])^2 / (rsq * sin_theta_sq)
        return cylinder.external_field * (2 * cos2 - 1) / rsq
    end
end


single_susceptibility_gradient(c::CylinderSusceptibility) = abs(c.external_field) / c.rsq^(3//2)


end
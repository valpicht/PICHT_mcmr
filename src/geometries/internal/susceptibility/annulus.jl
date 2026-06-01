"""
Off-resonance field due to magnetic susceptibility of an hollow cylinder.
"""
module Annulus

import StaticArrays: SVector
import ..Base: BaseSusceptibility, single_susceptibility, single_susceptibility_gradient

"""
    AnnulusSusceptibility(inner_radius, outer_radius, chi_I, chi_A, b0_field)

Creates an annular susceptibility source, which computes the off-resonance field expected for a myelin sheath.
"""
struct AnnulusSusceptibility <: BaseSusceptibility{2}
    inner_rsq :: Float64
    outer_rsq :: Float64
    internal_field :: Float64
    external_field :: Float64
    chi_I :: Float64
    chi_A :: Float64
end


function AnnulusSusceptibility(inner_radius::Number, outer_radius::Number, chi_I::Number, chi_A::Number, b0_field::AbstractVector)
    @assert length(b0_field) == 2
    sin_theta_sq = sum(b0_field .* b0_field)
    g_ratio = inner_radius / outer_radius
    AnnulusSusceptibility(
        inner_radius^2,
        outer_radius^2,
        -0.75 * chi_A * log(g_ratio) * sin_theta_sq,
        (chi_I + chi_A/4) * (outer_radius^2 - inner_radius^2) * sin_theta_sq / 2,
        chi_I,
        chi_A,
    )
end

"""
    single_susceptibility(annulus, position, distance, stuck_inside, b0_field)

Computed by the hollow cylinder fiber model from [Wharton_2012](@cite).
"""
function single_susceptibility(annulus::AnnulusSusceptibility, position::AbstractVector, distance::Number, stuck_inside::Union{Nothing, Bool}, b0_field::SVector{2, Float64})
    rsq = distance * distance

    if ~isnothing(stuck_inside) && rsq ≈ annulus.inner_rsq
        inside = stuck_inside ? 2 : 1
    elseif ~isnothing(stuck_inside) && rsq ≈ annulus.outer_rsq
        inside = stuck_inside ? 1 : 0
    else
        inside = rsq < annulus.inner_rsq ? 2 : (rsq < annulus.outer_rsq ? 1 : 0) 
    end

    if inside == 2
        return annulus.internal_field
    else
        sin_theta_sq = b0_field[1] * b0_field[1] + b0_field[2] * b0_field[2]
        if iszero(sin_theta_sq)
            cos2f = 0.5
        else
            cos2 = (b0_field[1] * position[1] + b0_field[2] * position[2])^2 / (rsq * sin_theta_sq)
            cos2f = 2 * cos2 - 1
        end
        if inside == 1
            return (
                annulus.chi_I * (2//3 - sin_theta_sq * (1 + cos2f * (annulus.inner_rsq / rsq))) / 2 +
                annulus.chi_A * (sin_theta_sq * (-5//12 - cos2f/8 * (1 + annulus.inner_rsq/rsq) + 3//8 * log(annulus.outer_rsq/rsq)) - (1-sin_theta_sq) / 6)
            )
        else
            return annulus.external_field * cos2f / rsq
        end
    end
end


single_susceptibility_gradient(a::AnnulusSusceptibility) = abs(a.external_field) / a.outer_rsq^(3//2)


end

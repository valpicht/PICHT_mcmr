"""
Base module for computing the off-resonance field from magnetic susceptibility sources.
"""
module Susceptibility
include("base.jl")
include("grid.jl")
include("cylinder.jl")
include("annulus.jl")
include("triangle.jl")

import .Base: BaseSusceptibility, single_susceptibility, single_susceptibility_gradient
import .Grid: FixedSusceptibility, SusceptibilityGrid, SusceptibilityGridElement, susceptibility_off_resonance, off_resonance_gradient, dipole_approximation, dipole_approximation_repeat, SusceptibilityGridNoRepeat, SusceptibilityGridRepeat, IsotropicSusceptibilityGridElement, AnisotropicSusceptibilityGridElement
import .Cylinder: CylinderSusceptibility
import .Annulus: AnnulusSusceptibility
import .Triangle: TriangleSusceptibility, IsotropicTriangleSusceptibility, AnisotropicTriangleSusceptibility, triangle_magnetisation

end
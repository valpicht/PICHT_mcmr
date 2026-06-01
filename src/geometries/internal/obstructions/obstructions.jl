"""
Defines all the base obstructions supported by MCMRSimulator.

They are all sub-types of `FixedObstruction`.
"""
module Obstructions
include("obstruction_intersections.jl")
include("fixed_obstructions.jl")
include("walls.jl")
include("rounds.jl")
include("OverlappingSphere.jl")
include("triangles.jl")
include("shifts.jl")

import .ObstructionIntersections: ObstructionIntersection, empty_obstruction_intersections, has_intersection
import .FixedObstructions: FixedObstruction, BoundingBox, isinside, has_inside, detect_intersection, obstruction_type, radius, size_scale, random_surface_positions
import .Walls: Wall
import .Rounds: Round, Cylinder, Sphere
import .Triangles: Triangle, FullTriangle, IndexTriangle, curvature
import .Shifts: Shift, get_quadrant
end
"""
Internal representation of the tissue geometry.

All objects are immutable and optimised for both memory usage and computational speed.
"""
module Internal
include("ray_grid_intersection.jl")
include("bounding_boxes.jl")
include("obstructions/obstructions.jl")
include("hit_grids.jl")
include("intersections.jl")
include("reflections.jl")
include("fixed_obstruction_groups.jl")
include("properties.jl")
include("susceptibility/susceptibility.jl")

import .FixedObstructionGroups: FixedObstructionGroup, FixedGeometry, repeating, isinside, detect_intersection, prepare_isinside!
import .Reflections: Reflection, direction, previous_hit, has_hit, has_intersection, empty_reflection
import .Intersections: Intersection, empty_intersection
import .HitGrids: HitGrid, detect_intersection_grid, obstructions

# Added `OverlappingSpheres` so the new internal geometry type is available
# through the `Internal` module. This is required for substrates represented
# as unions of overlapping spheres.
import .Obstructions: 
    FixedObstruction, ObstructionIntersection, empty_obstruction_intersections,
    Wall, Round, Cylinder, Sphere, OverlappingSpheres, Triangle, IndexTriangle, FullTriangle, Shift,
    has_inside, radius, obstruction_type, size_scale, get_quadrant, random_surface_positions

import .BoundingBoxes: BoundingBox, could_intersect, lower, upper
import .RayGridIntersection: ray_grid_intersections
import .Properties: R1, R2, off_resonance, permeability, surface_relaxation, surface_density, dwell_time, max_timestep_sticking, MRIProperties, max_permeability_non_inf, max_surface_relaxation, min_dwell_time
import .Susceptibility: FixedSusceptibility, susceptibility_off_resonance, off_resonance_gradient

end
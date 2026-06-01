"""
Defines the microstructural environment with which the spins interact.
"""
module Geometries
include("internal/internal.jl")
include("user/user.jl")
import .User: 
    fix, fix_susceptibility,
    ObstructionGroup, IndexedObstruction,
    Wall, Walls,
    Cylinder, Cylinders,
    Sphere, Spheres,
    Annulus, Annuli,
    Triangle, Mesh,
    Ring, BendyCylinder,
    load_mesh, random_positions_radii, nvolumes,
    write_geometry, read_geometry
end
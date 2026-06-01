"""
User interface for defining geometry.
"""
module User
include("obstructions/obstructions.jl")
include("split_mesh.jl")
include("size_scales.jl")
include("fix.jl")
include("fix_susceptibility.jl")
include("load_mesh.jl")
include("random_distribution.jl")
include("json.jl")
include("to_mesh.jl")
import .Fix: fix
import .FixSusceptibility: fix_susceptibility
import .Obstructions: ObstructionGroup, IndexedObstruction, nvolumes,
    Wall, Walls,
    Cylinder, Cylinders,
    Sphere, Spheres,
    Annulus, Annuli,
    Triangle, Mesh,
    Ring, BendyCylinder
import .SizeScales: size_scale
import .RandomDistribution: random_positions_radii
import .LoadMesh: load_mesh
import .JSON: write_geometry, read_geometry
end
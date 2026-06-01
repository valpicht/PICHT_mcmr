using MCMRSimulator
using MRIBuilder
pos, rad = random_positions_radii(30., 0.6, 2; variance=0.)
main_seq = DWI(bval=3., slice_thickness=2.)
geometries = Dict(
    "walls" => Walls(position=0.3, repeats=2.),
    "cylinders" => Cylinders(position=pos, radius=rad, repeats=(30., 30.)),
    "cylinders_myelin" => Cylinders(position=pos, radius=rad, repeats=(30., 30.), g_ratio=0.8),
    "cylinders_MT" => Cylinders(position=pos, radius=rad, repeats=(30., 30.), surface_density=1., dwell_time=30.),
    "cylinders_perm" => Cylinders(position=pos, radius=rad, repeats=(30., 30.), permeability=Inf),
    "bendy_cylinder" => BendyCylinder(control_point=[[0., 0., 0.], [0., 0.5, 1.], [0.5, 0., 2.]], radius=[0.6, 0.4, 0.5], repeats=[2., 2., 3.], closed=[0, 0, 1]),
)

mesh_fn = "/Users/michielcottaar/Work/data/microstructure-meshes/from_Palombo/cell_tissue/cell_tissue_replicated_kappa6.ply"
if isfile(mesh_fn)
    geometries["mesh"] = load_mesh(mesh_fn, repeats=(40, 40, 40))
    geometries["mesh_no_repeats"] = load_mesh(mesh_fn)
end

simulations = Dict{String, Simulation}(
    key => Simulation(main_seq; geometry=geometry, diffusivity=2., timestep=1e-2)
    for (key, geometry) in pairs(geometries)
)
simulations["no_diff"] = Simulation(main_seq; diffusivity=0., timestep=1e-2);
simulations["free_diff"] = Simulation(main_seq; diffusivity=0., timestep=1e-2);
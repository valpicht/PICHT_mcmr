# ============================================================
# utils_mesh.jl
#
# Utility functions for triangular mesh tests.
#
# This file is included by the scripts in:
#
#     docs/scripts/03_triangular_mesh_tests/
#
# It avoids hard-coded absolute paths and provides small helper
# functions reused across the mesh scripts.
# ============================================================

using MCMRSimulator
using MRIBuilder

"""
    default_mesh_dir()

Return the directory where example mesh files are stored.

The mesh files should be placed in:

    substrate_generation/example_outputs/

This avoids absolute paths such as:

    /home/valentine/project/Grow_neurons/results/...
"""
function default_mesh_dir()
    return joinpath(
        "substrate_generation",
        "example_outputs",
    )
end

"""
    add_ply_extension_if_missing(filename::String)

Add the `.ply` extension to a mesh name if it is missing.

Examples:

    add_ply_extension_if_missing("one_soma")
    # returns "one_soma.ply"

    add_ply_extension_if_missing("one_soma.ply")
    # returns "one_soma.ply"
"""
function add_ply_extension_if_missing(filename::String)
    if endswith(filename, ".ply")
        return filename
    else
        return filename * ".ply"
    end
end

"""
    mesh_path_from_name(mesh_name::String)

Return the relative path to a mesh file stored in the default mesh folder.

The input can be given with or without `.ply`.

Example:

    mesh_path_from_name("cell_sphere_r19_neurons_list_new_mesh")
"""
function mesh_path_from_name(mesh_name::String)
    filename = add_ply_extension_if_missing(mesh_name)
    return joinpath(default_mesh_dir(), filename)
end

"""
    load_example_mesh(mesh_name::String)

Load a triangular mesh from `substrate_generation/example_outputs/`.

Examples:

    geometry = load_example_mesh("cell_sphere_r19_neurons_list_new_mesh")
    geometry = load_example_mesh("cell_sphere_r19_neurons_list_new_mesh.ply")
"""
function load_example_mesh(mesh_name::String)
    mesh_path = mesh_path_from_name(mesh_name)

    println("Mesh path: ", mesh_path)
    @show isfile(mesh_path)

    if !isfile(mesh_path)
        error("""
        Mesh file not found.

        Expected file:
            $mesh_path

        Copy the mesh into:
            substrate_generation/example_outputs/
        """)
    end

    geometry = load_mesh(mesh_path)

    println("Loaded geometry type: ", typeof(geometry))

    return geometry
end

"""
    print_mesh_ranges(geometry)

Print the x, y, and z ranges of the first mesh object.
"""
function print_mesh_ranges(geometry)
    verts = geometry.vertices[1]

    println("Type of verts: ", typeof(verts))
    println("Number of vertices: ", length(verts))
    println("First vertex: ", verts[1])

    xmin = minimum(v[1] for v in verts)
    xmax = maximum(v[1] for v in verts)

    ymin = minimum(v[2] for v in verts)
    ymax = maximum(v[2] for v in verts)

    zmin = minimum(v[3] for v in verts)
    zmax = maximum(v[3] for v in verts)

    println("x range: ", xmin, " to ", xmax)
    println("y range: ", ymin, " to ", ymax)
    println("z range: ", zmin, " to ", zmax)
end

"""
    make_bbox_from_geometry(geometry; margin = 2.0)

Create a bounding box around the first mesh object.
"""
function make_bbox_from_geometry(geometry; margin = 2.0)
    verts = geometry.vertices[1]

    xmin = minimum(v[1] for v in verts)
    xmax = maximum(v[1] for v in verts)

    ymin = minimum(v[2] for v in verts)
    ymax = maximum(v[2] for v in verts)

    zmin = minimum(v[3] for v in verts)
    zmax = maximum(v[3] for v in verts)

    bbox = BoundingBox(
        [xmin - margin, ymin - margin, zmin - margin],
        [xmax + margin, ymax + margin, zmax + margin],
    )

    return bbox
end

"""
    create_default_dwi_sequence(; bval = 2.0, TE = 80, TR = 300)

Create a default DWI sequence used in the first mesh tests.
"""
function create_default_dwi_sequence(; bval = 2.0, TE = 80, TR = 300)
    return DWI(
        bval = bval,
        TE = TE,
        TR = TR,
        scanner = MRIBuilder.Siemens_Prisma,
    )
end
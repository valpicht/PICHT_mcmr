"""
Defines [`load_mesh`](@ref).
"""
module LoadMesh
import PlyIO
import ..Obstructions: Mesh

"""
    ply_from_mesh(file)

Loads a [`Mesh`](@ref) from a PLY file.

PLY stands for Polygon File Format (http://paulbourke.net/dataformats/ply/).
PLY IO is handled by PlyIO.jl (https://github.com/JuliaGeometry/PlyIO.jl).
"""
function ply_from_mesh(ply_file; kwargs...)
    ply = PlyIO.load_ply(ply_file)
    vertices = collect(zip(
        Float32.(ply["vertex"]["x"]),
        Float32.(ply["vertex"]["y"]),
        Float32.(ply["vertex"]["z"]),
    ))
    indices = [v .+ 1 for v in ply["face"]["vertex_indices"]]
    return Mesh(vertices=vertices, triangles=indices; kwargs...)
end


"""
    load_mesh(file)

Loads a [`Mesh`](@ref) from a file.

Currently only PLY files are supported.

## PLY format
PLY stands for Polygon File Format (http://paulbourke.net/dataformats/ply/).
PLY IO is handled by PlyIO.jl (https://github.com/JuliaGeometry/PlyIO.jl).
"""
function load_mesh(io::IO; kwargs...)
    mark(io)
    firstline = readline(io)
    reset(io)
    if firstline == "ply"
        return ply_from_mesh(io; kwargs...)
    else
        throw(ErrorException("Could not recognise input mesh format!"))
    end
end


function load_mesh(file_name::AbstractString; kwargs...)
    open(file_name, "r") do fid
        load_mesh(fid; kwargs...)
    end
end

end
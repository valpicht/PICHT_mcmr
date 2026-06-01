"""
Mesh-specific operations to normalise (i.e., ensure normals point outwards) and split a mesh.
"""
module SplitMesh

import StaticArrays: SVector, MVector
import SparseArrays: sparse, SparseMatrixCSC
import LinearAlgebra: norm, ⋅
import Statistics: mean
import ...Internal.Obstructions.Triangles: normal, curvature
import ..Obstructions: Mesh, isglobal, nvolumes


"""
    fix_mesh(mesh)

Returns a fixed mesh, where all normals point outwards.
"""
function fix_mesh(old_mesh::Mesh)
    new_triangles = MVector{3, Int}.(old_mesh.triangles.value)
    triangle_indices = components(old_mesh)
    for index in unique(components(old_mesh))
        triangles = new_triangles[triangle_indices .== index]
        make_normals_consistent!(triangles)
        if curvature(triangles, old_mesh.vertices.value) < 0
            # flip all triangles
            for t in triangles
                v = t[1]
                t[1] = t[2]
                t[2] = v
            end
        end
    end
    kwargs = Dict{Symbol, Any}(key => getproperty(old_mesh, key).value for key in old_mesh.unique_keys)
    return Mesh(; number=length(new_triangles), kwargs..., triangles=new_triangles)
end

"""
    make_normals_consistent!(triangles)

Adjust the triangles to all point outwards or all point inwards.
Assumes that all the triangles are connected (can be enforced using [`connected_components`](@ref)).
"""
function make_normals_consistent!(triangles::AbstractVector)
    edges(t) = ((t[1], t[2]), (t[2], t[3]), (t[3], t[1]))
    counter = Dict{Tuple{Int, Int}, Int}()
    for triangle in triangles
        for (index1, index2) in edges(triangle)
            edge = index2 > index1 ? (index1, index2) : (index2, index1)
            if edge in keys(counter)
                counter[edge] += 1
            else
                counter[edge] = 1
            end
        end
    end
    triangles_seen = fill(false, length(triangles))
    edges_seen = Set{Tuple{Int, Int}}()
    while ~all(triangles_seen)
        starting_triangle = findfirst(.~triangles_seen)
        triangles_seen[starting_triangle] = true
        nseen = 0
        while sum(triangles_seen) > nseen
            nseen = sum(triangles_seen)
            union!(edges_seen, edges(triangles[starting_triangle]))
            for (i, triangle) in enumerate(triangles)
                if triangles_seen[i]
                    continue
                end
                flip_required = -1
                for (index1, index2) in edges(triangle)
                    norm_edge = index2 > index1 ? (index1, index2) : (index2, index1)
                    if counter[norm_edge] != 2
                        continue
                    end
                    if (index1, index2) in edges_seen
                        @assert flip_required in (-1, 1)
                        flip_required = 1
                    elseif (index2, index1) in edges_seen
                        @assert flip_required in (-1, 0)
                        flip_required = 0
                    end
                end
                if flip_required == -1
                    continue
                elseif flip_required == 1
                    v = triangle[1]
                    triangle[1] = triangle[2]
                    triangle[2] = v
                end
                triangles_seen[i] = true
                union!(edges_seen, edges(triangle))
            end
        end
    end
end

"""
    connectivity_matrix(triangles, nvertices)

Create a nverticesxnvertices sparse boolean matrix, which is true for any connected vertices.
"""
function connectivity_matrix(triangles)
    t1 = [t[1] for t in triangles]
    t2 = [t[2] for t in triangles]
    t3 = [t[3] for t in triangles]
    nvertices = max(maximum(t1), maximum(t2), maximum(t3))
    sparse(vcat(t1, t2, t3), vcat(t2, t3, t1), true, nvertices, nvertices)
end


"""
    connected_indices(connectivity_matrix)
    connected_indices(triangles)

Returns a vector of indices, where each connected component has been given a unique index (starting with 1).

For the triangles, each vertex will be given an index (not each triangle).
"""
connected_indices(t::AbstractVector) = connected_indices(connectivity_matrix(t))

function connected_indices(m::SparseMatrixCSC)
    nvertices = size(m, 1)
    not_assigned = fill(true, nvertices)
    result = fill(0, nvertices)
    component_index = 1
    while any(not_assigned)
        start_index :: Int = findfirst(not_assigned)
        new_component = BitVector(fill(false, nvertices))
        tocheck = BitVector(fill(false, nvertices))

        new_component[start_index] = true
        tocheck[start_index] = true

        repeat = true
        while repeat
            repeat = false
            for i in 1:nvertices
                if @inbounds tocheck[i]
                    for j in m.rowval[m.colptr[i]:m.colptr[i + 1] - 1]
                        if ~new_component[j]
                            new_component[j] = true
                            tocheck[j] = true
                            if j < i
                                repeat = true
                            end
                        end
                    end
                    tocheck[i] = false
                end
            end
        end

        @assert sum(not_assigned .& new_component) == sum(new_component) "new component includes elements already assigned to previous components"
        not_assigned[new_component] .= false
        result[new_component] .= component_index
        component_index += 1
    end
    @assert all(result .> 0)
    return result
end


"""
    components(mesh)

Returns which component each element in the mesh belongs to.

If not set explicitly (using `mesh.components=[...]`), it will be computed once based on the connectivity structure.
"""
function components(mesh::Mesh)
    if mesh.n_obstructions == 0
        return Int[]
    end
    if isnothing(mesh.components.value)
        vertex_indices = connected_indices(mesh.triangles.value)
        mesh.components.value = [vertex_indices[t[1]] for t in mesh.triangles.value]
    end
    return mesh.components.value
end

nvolumes(mesh::Mesh) = maximum(components(mesh); init=0)

end
"""
Defines `ray_grid_intersections`, which computes the intersections between a spin path and some grid.
"""
module RayGridIntersection

import StaticArrays: SVector

struct RayGridIntersections{N}
    origin :: SVector{N, Float64}
    destination :: SVector{N, Float64}
    direction :: SVector{N, Float64}
    abs_inv_direction :: SVector{N, Float64}
end

"""
    ray_grid_intersections(origin, destination)

Computes all voxels crossed by a ray between `origin` and `destination` with an inifinitely extending 1x1x1 grid.
Both origin and destination are length-3 vectors.
The returned object is an iterator returning a tuple with:
- N-length vector with the voxel that we are crossing through
- Float with the time the ray entered voxel (0=`origin`, 1=`destination`)
- N-length vector with position within voxel that the ray entered (i.e., numbers between 0 and 1)
- Float with the time the ray left the voxel (0=`origin`, 1=`destination`)
- N-length vector with position within voxel that the ray left (i.e., numbers between 0 and 1)
"""
function ray_grid_intersections(origin :: SVector{N, Float64}, destination :: SVector{N, Float64}) where {N}
    direction = destination - origin
    RayGridIntersections(origin, destination, direction, map(d -> Float64(1/abs(d)), direction))
end

function ray_grid_intersections(origin :: AbstractVector, destination :: AbstractVector)
    N = length(origin)
    ray_grid_intersections(SVector{N, Float64}(origin), SVector{N, Float64}(destination))
end

Base.iterate(rgi::RayGridIntersections) = Base.iterate(rgi, (rgi.origin, zero(Float64), map(o -> Int(floor(o)), rgi.origin)))
function Base.iterate(rgi::RayGridIntersections{N}, state::Tuple{SVector{N, Float64}, Float64, SVector{N, Int}}) where {N}
    (prev_pos::SVector{N, Float64}, prev_time::Float64, current_voxel::SVector{N, Int}) = state
    if prev_time >= 1.
        return nothing
    end
    time_to_hit::Float64, dimension::Int = next_hit(prev_pos, current_voxel, rgi.direction, rgi.abs_inv_direction)
    next_time = prev_time + time_to_hit
    if next_time > 1.
        return (
            (current_voxel, prev_time, prev_pos - current_voxel, one(Float64), rgi.destination - current_voxel), 
            (rgi.destination, next_time, current_voxel)
        )
    end
    ddim = rgi.direction[dimension] > 0 ? 1 : -1
    next_voxel = map((i, v) -> i == dimension ? v + ddim : v, 1:N, current_voxel)
    next_pos = rgi.origin .+ (rgi.direction .* next_time)
    res = (
        (current_voxel, prev_time, prev_pos - current_voxel, next_time, next_pos - current_voxel),
        (next_pos, next_time, next_voxel)
    )
    return res
end

function next_hit(prev_pos::SVector{N, Float64}, current_voxel::SVector{N, Int}, direction::SVector{N, Float64}, abs_inv_direction::SVector{N, Float64}) where {N}
    time_to_hit = Float64(Inf)
    dimension = 0
    for dim in 1:N
        within_voxel = prev_pos[dim] - current_voxel[dim]
        d = direction[dim]
        if iszero(d)
            continue
        end
        next_hit = (d > 0. ? 1. - within_voxel : within_voxel) * abs_inv_direction[dim]
        if next_hit < time_to_hit
            time_to_hit = next_hit
            dimension = dim
        end
    end
    return time_to_hit, dimension
end

Base.length(rgi::RayGridIntersections) = sum(abs.(Int.(floor.(rgi.destination)) .- Int.(floor.(rgi.origin)))) + 1
Base.eltype(::RayGridIntersections{N}) where {N} = Tuple{SVector{N, Float64}, Float64, SVector{N, Float64}, Float64, SVector{N, Float64}}

end
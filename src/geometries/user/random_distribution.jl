"""
Defines [`random_positions_radii`](@ref).
"""
module RandomDistribution
import StaticArrays: SVector, MVector
import Distributions

function draw_radii(box_size::SVector{N, <:Real}, target_density::Real, distribution::Distributions.ContinuousUnivariateDistribution, allowed_range::Tuple{<:Real, <:Real}) where {N}
    if target_density > 1 || target_density < 0
        error("Target density should be a number between 0 and 1")
    end
    @assert allowed_range[2] > allowed_range[1]
    volume = 0.
    target_volume = prod(box_size) * target_density
    samples = Float64[]
    while volume < target_volume
        new_value = rand(distribution)
        if new_value < allowed_range[1] || new_value > allowed_range[2]
            continue 
        end
        push!(samples, new_value)
        if N == 3
            volume += 4//3 * π * new_value^3
        elseif N == 2
            volume += π * new_value^2
        elseif N == 1
            volume += 2 * new_value
        end
    end
    samples
end


norm_dist(orig_dist::Float64, half_repeat::Float64) = mod(orig_dist + half_repeat, half_repeat * 2) - half_repeat

function eliminate_overlap!(positions::AbstractVector{SVector{N, Float64}}, radii::AbstractVector{Float64}, box_size::SVector{N, Float64}; max_iter=1000) where {N}
    @assert length(positions) == length(radii)
    moved = true
    half_box_size = box_size ./ 2

    step_size = 1.

    for _ in 1:max_iter
        moved = false
        for i in eachindex(positions)
            for j in 1:(i-1)
                displacement = map(norm_dist, positions[i] - positions[j], half_box_size)
                distsq = sum(p->p*p, displacement)
                min_dist = radii[i] + radii[j]
                if distsq < (min_dist * min_dist)
                    moved = true
                    dist = sqrt(distsq)
                    total_movement = max((min_dist - dist) * (1 + step_size), 0.01 * min_dist)
                    # move larger cylinder less
                    positions[i] = positions[i] .+ displacement .* (total_movement * radii[j] / (min_dist * dist))
                    positions[j] = positions[j] .- displacement .* (total_movement * radii[i] / (min_dist * dist))
                end
            end
        end
        step_size *= 0.8
        if !moved
            positions[:] = map(pos -> norm_dist.(pos, half_box_size), positions)
            return
        end
    end
    error("No solution without overlap found after max_iter=$max_iter iterations")
end


function repel_position!(positions::AbstractVector{SVector{N, Float64}}, radii::AbstractVector{Float64}, box_size::SVector{N, Float64}; max_iter=1000, repulsion_strength=1e-3) where {N}

    half_box_size = box_size ./ 2

    function does_repulsed_have_a_collision(proposed_pos, index, positions)

        for i in eachindex(positions)
            if i != index
                min_dist_sq = (radii[index] + radii[i])^2
                dist_sq = sum(norm_dist.(proposed_pos .- positions[i], half_box_size) .^ 2)
                if min_dist_sq >= dist_sq
                    return true
                end
            end
        end
        return false
    end
    
    proposed_pos = zero(MVector{N, Float64})
    # Update circle positions based on repulsion
    for _ in 1:max_iter
        success_count=0
        for i in eachindex(positions)
            proposed_pos[:] = positions[i]
            for j in eachindex(positions)
                if i != j
                    displacement = norm_dist.(positions[j] .- positions[i], half_box_size)

                    dist_sq = sum(displacement .* displacement)

                    # force = (radius^2*radius^2)*repulsion_strength / (dist * mean_radius^2)
                    force = (radii[i]^N * radii[j]^N)*repulsion_strength / (dist_sq^(N//2) * (radii[i] + radii[j])^N)
                    proposed_pos .-= force .* displacement
                end
            end
            # Abandon changes if there's overlap
            if !does_repulsed_have_a_collision(proposed_pos, i, positions)
                positions[i] = proposed_pos
                success_count += 1
            end
        end
        if iszero(success_count)
            # Stop early when equilibrium is reached and nothing can be moved in a loop
            break
        end
    end
end


"""
    random_positions_radii(box_size, target_density, n_dimensions; distribution=Gamma, mean=1., variance=1., max_iter=1000, min_radius=0.1, max_radius=Inf)

Randomly distributes circles or spheres in space.

Arguments:
- `box_size`: Size of the infinitely repeating box of random positions
- `target_density`: Final density of the circles/spheres. This density will only be approximately reached
- `n_dimensions`: dimensionality of the space (2 for cicles; 3 for spheres)
- `distribution`: distribution from which the radii are drawn (from [Distributions.jl](https://juliastats.org/Distributions.jl/stable/))
- `mean`: mean of the gamma distribution (ignored if `distribution` explicitly set)
- `variance`: variance of the gamma distribution (ignored if `distribution` explicitly set)
- `max_iter`: maximum number of iterations to try to prevent the circles/spheres from overlapping. An error is raised if they still overlap after this number of iterations.
- `repulsion_strength`: strength of the repulsion in each step (defaults to 0.001).
- `max_iter_repulse`: maximum number of iterations that circles/spheres will be repulsed from each other
- `min_radius`: samples from the distribution below this radius will be rejected (in um).
- `max_radius`: samples from the distribution above this radius will be rejected (in um).
"""
function random_positions_radii(
    box_size, target_density::Real, ndim::Int;
    distribution=nothing, mean=1., variance=0.1, max_iter=1000,
    repulsion_strength=1e-3, max_iter_repulse=1000,
    min_radius=0.1, max_radius=Inf,
)
    if iszero(variance)
        distribution = Distributions.Normal(mean, 0.)
    elseif isnothing(distribution)
        a = mean * mean / variance
        scale = variance/mean
        distribution = Distributions.Gamma(a, scale)
    end
    if isa(box_size, Real)
        box_size = SVector{ndim, Float64}(fill(box_size, ndim))
    else
        box_size = SVector{ndim, Float64}(box_size)
    end
    target_density = Float64(target_density)
    radii = draw_radii(box_size, target_density, distribution, (min_radius, max_radius))
    positions = [SVector{ndim, Float64}((rand(ndim) .- 0.5) .* box_size) for _ in 1:length(radii)]
    eliminate_overlap!(positions, radii, box_size; max_iter=max_iter)
    if !isone(ndim)
        repel_position!(positions, radii, box_size; max_iter=max_iter_repulse, repulsion_strength=repulsion_strength)
    end
    return positions, radii
end

end
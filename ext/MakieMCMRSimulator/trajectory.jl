module Trajectory
using Makie
import Colors
import MCMRSimulator.Plot: PlotPlane, project_trajectory, Plot_Trajectory
import ..Utils: Utils
import MCMRSimulator.Spins: Snapshot, position, get_sequence


function Makie.plot!(scene::Plot_Trajectory{<:Tuple{Nothing, <:AbstractVector{<:Snapshot}}})
    Makie.@extract scene (trajectory, color, sequence)
    nspins = @lift length($trajectory[1])
    kwargs = Dict([key => scene[key] for key in [
        :visible, :overdraw, :fxaa, :transparency, :inspectable, :depth_shift, :model, :space,
    ]])
    lift(nspins) do N 
        for index in 1:N
            positions = @lift map(s -> Makie.Point3f(position(s[index])), $trajectory)
            colors = @lift $color === Makie.automatic ? map(s -> Utils.color(s[index]; sequence=$sequence), $trajectory) : $color
            Makie.lines!(scene, positions; color=colors, kwargs...)
        end
    end
    scene
end

Makie.plottype(::AbstractVector{<:Snapshot}) = Plot_Trajectory


function Makie.plot!(scene::Plot_Trajectory{<:Tuple{<:PlotPlane, <:AbstractVector{<:Snapshot}}})
    function _get_colors(sequence_index :: Integer, spin_index :: Integer, snapshots :: Vector{<:Snapshot{N}}, times :: Vector{Float64}) where {N}
        if N == 0 || sequence_index == 0
            return [Colors.HSV() for _ in times]
        else
            colors_main = map(s -> Utils.color(get_sequence(s[spin_index], sequence_index)), snapshots)
            return [isfinite(time) ? colors_main[Int(round(time))] : Colors.HSV() for time in times]
        end
    end
    plot_plane = scene[1]
    trajectory = scene[2]
    Makie.@extract scene (color, sequence)
    kwargs = Dict([key => scene[key] for key in [
        :visible, :overdraw, :fxaa, :transparency, :inspectable, :depth_shift, :model, :space,
    ]])
    nspins = @lift length($trajectory[1])
    lift(nspins) do N 
        for index in 1:N
            positions_3d = @lift map(s -> position(s[index]), $trajectory)
            projected = @lift project_trajectory($plot_plane, $positions_3d)
            positions = @lift Makie.Point2f.($projected[1])
            colors = @lift $color == Makie.automatic ? _get_colors($(scene[:sequence]), index, $trajectory, $projected[2]) : $color
            Makie.lines!(scene, positions; color=colors, kwargs...)
        end
    end
    scene
end

Makie.plottype(::PlotPlane, ::AbstractVector{<:Snapshot}) = Plot_Trajectory

end
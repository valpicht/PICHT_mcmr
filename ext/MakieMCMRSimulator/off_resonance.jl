module OffResonance
using Makie
import StaticArrays: SVector
import MCMRSimulator.Plot: PlotPlane, Plot_Off_Resonance
import MCMRSimulator.Geometries.Internal: FixedSusceptibility, susceptibility_off_resonance
import MCMRSimulator.Geometries: fix_susceptibility

function Makie.plot!(scene::Plot_Off_Resonance)
    plot_plane = scene[1]
    raw_geometry = scene[2]

    kwargs = Dict([key => scene[key] for key in [
        :visible, :overdraw, :fxaa, :transparency, :inspectable, :depth_shift, :model, :space,
        :colormap, :colorscale, :colorrange, :nan_color, :lowclip, :highclip, :alpha,
    ]])
    susc = @lift raw_geometry isa FixedSusceptibility ? raw_geometry : fix_susceptibility($raw_geometry)

    dims = @lift -0.5:(1/$(scene[:ngrid])):0.5
    xx_1d = @lift $dims * $plot_plane.sizex
    yy_1d = @lift $dims * $plot_plane.sizey
    pos_plane = @lift broadcast(
        (x, y) -> SVector{3}([x, y, 0.]),
        reshape($xx_1d, length($xx_1d), 1),
        reshape($yy_1d, 1, length($yy_1d)),
    )
    pos_orig = @lift inv($plot_plane.transformation).($pos_plane)
    field = @lift map(p->susceptibility_off_resonance($susc, p), $pos_orig)
    x_interval = @lift (-0.5 * $plot_plane.sizex) .. (0.5 * $plot_plane.sizex)
    y_interval = @lift (-0.5 * $plot_plane.sizey) .. (0.5 * $plot_plane.sizey)
    Makie.image!(scene, x_interval, y_interval, field; kwargs...)
    scene
end

end
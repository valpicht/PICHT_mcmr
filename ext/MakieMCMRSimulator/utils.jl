module Utils
import Makie
import Colors
import MCMRSimulator.Spins: Spin, SpinOrientation, phase, transverse
import MCMRSimulator: Snapshot
import MCMRSimulator.Plot: PlotPlane, print_sequence, Projectable
import MRIBuilder: plot_sequence


"""
    color(orient::SpinOrientation; saturation=1.)
    color(spin::Spin; saturation=1., sequence=1)

Returns a color representing the spin orientation in the transverse (x-y) plane.
Brighter colors have a larger transverse component, so that spins with no transverse component are black.
The actual color encodes the spin orientation.
"""
color(orient::SpinOrientation; saturation=1.) = Colors.HSV(phase(orient) + 180, saturation, transverse(orient))
color(spin::Spin{N}; sequence=1, kwargs...) where {N} = color(spin.orientations[sequence]; kwargs...)


Makie.args_preferred_axis(::Projectable) = Makie.LScene
Makie.args_preferred_axis(::PlotPlane) = Makie.Axis

function print_sequence(; sequence_file, output_file, t0, t1, kwargs...)
    sequence = read_sequence(sequence_file)
    f = plot_sequence(sequence)
    xlims!(f.axis, t0, t1)
    Makie.save(output_file, f)
end
end
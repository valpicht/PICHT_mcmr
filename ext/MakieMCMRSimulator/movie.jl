module Movie
using Makie
import MCMRSimulator.Plot: PlotPlane, simulator_movie
import MCMRSimulator.Spins: transverse
import MCMRSimulator.Simulations: Simulation
import MCMRSimulator.Evolve: readout

function simulator_movie(filename, simulator::Simulation{N}, times, repeats; figure_size=(1600, 800), trajectory_init=30, signal_init=10000, framerate=50, plane_orientation=:z, lengthscale=1., arrowsize=10.) where {N}
    if isa(trajectory_init, Integer)
        trajectory_init = [rand(3) .* repeats .- repeats ./ 2 for _ in 1:trajectory_init]
    end
    sig = readout(signal_init, simulator, times)
    if ndims(sig) == 1
        trans = [transverse.(sig) ./ signal_init]
    else
        trans = [transverse.(sig[index, :]) ./ signal_init for index in 1:N]
    end
    traj = readout(trajectory_init, simulator, times, return_snapshot=true);
    if ndims(traj) > 1
        traj = traj[1, :]
    end
    pp = PlotPlane(plane_orientation, sizex=repeats[1], sizey=repeats[2])

    index = Observable(1)
    time = @lift times[$index]
    fig = Figure(size=figure_size)
    ax_walk = Axis(fig[1, 1], title=(@lift("time = $(round($time, digits=1)) ms")))
    plot!(ax_walk, pp, simulator.geometry)
    xlims!(ax_walk, -repeats[1]/2, repeats[1]/2)
    ylims!(ax_walk, -repeats[2]/2, repeats[2]/2)
    plot!(ax_walk, pp, @lift(traj[$index]), lengthscale=lengthscale, arrowsize=arrowsize, kind=:dyad)

    ax_seq = Axis(fig[1, 2])

    function set_nan(arr, index)
        new = copy(arr)
        new[index+1:end] .= NaN
        new
    end

    plot!(ax_seq, simulator.sequences[1])
    for t in trans
        lines!(ax_seq, times, (@lift(set_nan(t, $index))))
    end
    xlims!(ax_seq, -3, maximum(times) + 5)


    record(fig, filename, 1:length(times);
            framerate=framerate) do i
        index[] = i
    end
end

end
# ============================================================
# 05_signal_vs_b_TE_and_ADC.jl
#
# Goal:
# Simulate triangular mesh signals over several b-values and echo
# times, compute ln(S/S0), fit a simple ADC model, and save plots.
#
# This script combines and cleans the logic from:
#   - run_simulation.jl
#   - plot_ln_signal_vs_b_TE.jl
#   - plot_adc_from_csv.jl
#
# Model:
#
#     ln(S/S0) = c - ADC * b
#
# Run from the repository root with:
#
#     julia +1.11 --project=. docs/scripts/03_triangular_mesh_tests/05_signal_vs_b_TE_and_ADC.jl
#
# Or run cell by cell in VS Code.
# ============================================================


# %% Cell 1 — Load packages and utilities

using MCMRSimulator
using MRIBuilder
using CairoMakie
using DelimitedFiles
using LsqFit

include(joinpath(@__DIR__, "utils_mesh.jl"))


# %% Cell 2 — Choose mesh and output paths

const MESH_NAME = "one_soma_list_new_mesh"

const OUTPUT_DIR = joinpath("docs", "data", "processed")
const FIGURE_DIR = joinpath("docs", "figures")

mkpath(OUTPUT_DIR)
mkpath(FIGURE_DIR)

const SIGNAL_CSV = joinpath(OUTPUT_DIR, "mesh_ln_signal_vs_b_TE_data.csv")
const ADC_CSV = joinpath(OUTPUT_DIR, "mesh_adc_vs_TE_data.csv")

const INSIDE_FIT_FIG = joinpath(FIGURE_DIR, "mesh_inside_ln_signal_vs_b_TE_fit.png")
const OUTSIDE_FIT_FIG = joinpath(FIGURE_DIR, "mesh_outside_ln_signal_vs_b_TE_fit.png")
const ADC_FIG = joinpath(FIGURE_DIR, "mesh_adc_vs_TE.png")

println("Selected mesh: ", MESH_NAME)


# %% Cell 3 — Define simulation helper

"""
    simulate_signal_vs_b(spins, geometry, bvals, TE; TR = 300, diffusivity = 1.0)

Run one DWI simulation per b-value for one echo time TE.

Returns:
- a vector of scalar transverse signal values.
"""
function simulate_signal_vs_b(
    spins,
    geometry,
    bvals,
    TE;
    TR = 300,
    diffusivity = 1.0,
)
    signals = Float64[]

    for b in bvals
        println("Running simulation for TE = $TE ms, b = $b")

        try
            sequence = DWI(
                bval = b,
                TE = TE,
                TR = TR,
                scanner = MRIBuilder.Siemens_Prisma,
            )

            sim = Simulation(
                sequence,
                diffusivity = diffusivity,
                geometry = geometry,
            )

            sig = readout(
                spins,
                sim,
                skip_TR = 2,
            )

            # Original scripts used sig.orient.transverse.
            # abs(...) makes the scalar robust if the transverse signal is complex.
            push!(signals, abs(sig.orient.transverse))

        catch err
            println("WARNING: failed for TE = $TE ms, b = $b")
            println("Error: ", err)
            push!(signals, NaN)
        end
    end

    return signals
end


# %% Cell 4 — Define log-normalization helper

"""
    log_normalized_signal(signals)

Compute ln(S/S0) for a vector of signal values.

If a signal is NaN or non-positive, the function returns NaN for that
point. This avoids crashing the ADC fit.
"""
function log_normalized_signal(signals)
    s0 = signals[1]
    ln_norm = Float64[]

    for s in signals
        if !isnan(s) && !isnan(s0) && s > 0 && s0 > 0
            push!(ln_norm, log(s / s0))
        else
            push!(ln_norm, NaN)
        end
    end

    return ln_norm
end


# %% Cell 5 — Define ADC model and fitting helper

"""
    adc_model(b, p)

Linear ADC model:

    ln(S/S0) = c - ADC * b

where:
- p[1] = c
- p[2] = ADC
"""
function adc_model(b, p)
    c, ADC = p
    return c .- ADC .* b
end

"""
    fit_adc_curvefit(bvals, ln_values)

Fit the ADC model to log-normalized signal values.

Returns:
- ADC estimate;
- fitted intercept;
- fitted y-values for plotting.
"""
function fit_adc_curvefit(bvals, ln_values)
    valid = .!isnan.(bvals) .& .!isnan.(ln_values)

    bfit = bvals[valid]
    yfit = ln_values[valid]

    if length(bfit) < 2
        return NaN, NaN, fill(NaN, length(bvals))
    end

    p0 = [0.0, 1.0]   # initial guess: intercept ≈ 0, ADC ≈ 1
    fit = curve_fit(adc_model, bfit, yfit, p0)

    c_fit = fit.param[1]
    adc_fit = fit.param[2]

    y_pred = adc_model(bvals, fit.param)

    return adc_fit, c_fit, y_pred
end


# %% Cell 6 — Load mesh

geometry = load_example_mesh(MESH_NAME)

println("Mesh loaded.")


# %% Cell 7 — Create bounding box and split spins

bbox = make_bbox_from_geometry(
    geometry;
    margin = 2.0,
)

println("Bounding box:")
println(bbox)

all_spins = Snapshot(10_000, bbox)

inside_mask = isinside(geometry, all_spins) .> 0
inside_spins = all_spins[inside_mask]
outside_spins = all_spins[.!inside_mask]

println("Number of inside spins: ", length(inside_spins))
println("Number of outside spins: ", length(outside_spins))


# %% Cell 8 — Define b-values and TE values

bvals = collect(0.0:0.1:1.0)

# Realistic TE values used in the original script.
TE_vals = [20, 40, 80]

println("b-values:")
println(bvals)

println("TE values:")
println(TE_vals)


# %% Cell 9 — Prepare containers

inside_ln_dict = Dict{Int, Vector{Float64}}()
outside_ln_dict = Dict{Int, Vector{Float64}}()

inside_fit_dict = Dict{Int, Vector{Float64}}()
outside_fit_dict = Dict{Int, Vector{Float64}}()

inside_adc = Float64[]
outside_adc = Float64[]

inside_intercepts = Float64[]
outside_intercepts = Float64[]

csv_rows = Any[]
push!(
    csv_rows,
    [
        "TE",
        "b",
        "inside_signal",
        "inside_ln_S_over_S0",
        "outside_signal",
        "outside_ln_S_over_S0",
    ],
)


# %% Cell 10 — Run simulations and fit ADC

for TE in TE_vals
    println("\n==============================")
    println("Processing TE = $TE ms")
    println("==============================")

    inside_signals = simulate_signal_vs_b(
        inside_spins,
        geometry,
        bvals,
        TE;
        TR = 300,
        diffusivity = 1.0,
    )

    outside_signals = simulate_signal_vs_b(
        outside_spins,
        geometry,
        bvals,
        TE;
        TR = 300,
        diffusivity = 1.0,
    )

    inside_ln = log_normalized_signal(inside_signals)
    outside_ln = log_normalized_signal(outside_signals)

    inside_ln_dict[TE] = inside_ln
    outside_ln_dict[TE] = outside_ln

    inside_adc_fit, inside_c_fit, inside_y_pred = fit_adc_curvefit(
        bvals,
        inside_ln,
    )

    outside_adc_fit, outside_c_fit, outside_y_pred = fit_adc_curvefit(
        bvals,
        outside_ln,
    )

    inside_fit_dict[TE] = inside_y_pred
    outside_fit_dict[TE] = outside_y_pred

    push!(inside_adc, inside_adc_fit)
    push!(outside_adc, outside_adc_fit)

    push!(inside_intercepts, inside_c_fit)
    push!(outside_intercepts, outside_c_fit)

    println("Inside fit for TE = $TE ms:")
    println("   ln(S/S0) = $inside_c_fit - $inside_adc_fit * b")

    println("Outside fit for TE = $TE ms:")
    println("   ln(S/S0) = $outside_c_fit - $outside_adc_fit * b")

    for i in eachindex(bvals)
        push!(
            csv_rows,
            [
                TE,
                bvals[i],
                inside_signals[i],
                inside_ln[i],
                outside_signals[i],
                outside_ln[i],
            ],
        )
    end
end


# %% Cell 11 — Save signal CSV

open(SIGNAL_CSV, "w") do io
    writedlm(io, csv_rows, ',')
end

println("\nSignal data saved:")
println(" - ", SIGNAL_CSV)


# %% Cell 12 — Save ADC summary CSV

adc_summary = hcat(
    TE_vals,
    inside_adc,
    inside_intercepts,
    outside_adc,
    outside_intercepts,
)

open(ADC_CSV, "w") do io
    writedlm(
        io,
        ["TE" "inside_ADC" "inside_intercept" "outside_ADC" "outside_intercept"],
        ',',
    )
    writedlm(io, adc_summary, ',')
end

println("ADC summary saved:")
println(" - ", ADC_CSV)


# %% Cell 13 — Plot inside ln(S/S0) with fitted curves

fig_inside = Figure(size = (900, 600))
ax_inside = Axis(
    fig_inside[1, 1],
    title = "Inside mesh: ln(S(b)/S(0)) vs b",
    xlabel = "b",
    ylabel = "ln(S(b)/S(0))",
)

for TE in TE_vals
    scatter!(
        ax_inside,
        bvals,
        inside_ln_dict[TE],
        label = "data TE = $TE ms",
    )

    lines!(
        ax_inside,
        bvals,
        inside_fit_dict[TE],
        label = "fit TE = $TE ms",
    )
end

axislegend(ax_inside, position = :rb)

save(INSIDE_FIT_FIG, fig_inside)
display(fig_inside)

println("Inside fit figure saved:")
println(" - ", INSIDE_FIT_FIG)


# %% Cell 14 — Plot outside ln(S/S0) with fitted curves

fig_outside = Figure(size = (900, 600))
ax_outside = Axis(
    fig_outside[1, 1],
    title = "Outside mesh: ln(S(b)/S(0)) vs b",
    xlabel = "b",
    ylabel = "ln(S(b)/S(0))",
)

for TE in TE_vals
    scatter!(
        ax_outside,
        bvals,
        outside_ln_dict[TE],
        label = "data TE = $TE ms",
    )

    lines!(
        ax_outside,
        bvals,
        outside_fit_dict[TE],
        label = "fit TE = $TE ms",
    )
end

axislegend(ax_outside, position = :rb)

save(OUTSIDE_FIT_FIG, fig_outside)
display(fig_outside)

println("Outside fit figure saved:")
println(" - ", OUTSIDE_FIT_FIG)


# %% Cell 15 — Plot ADC vs TE

fig_adc = Figure(size = (900, 600))
ax_adc = Axis(
    fig_adc[1, 1],
    title = "ADC vs TE",
    xlabel = "TE (ms)",
    ylabel = "ADC",
)

lines!(
    ax_adc,
    TE_vals,
    inside_adc,
    label = "Inside ADC",
)
scatter!(ax_adc, TE_vals, inside_adc)

lines!(
    ax_adc,
    TE_vals,
    outside_adc,
    label = "Outside ADC",
)
scatter!(ax_adc, TE_vals, outside_adc)

axislegend(ax_adc, position = :rb)

save(ADC_FIG, fig_adc)
display(fig_adc)

println("ADC figure saved:")
println(" - ", ADC_FIG)


# %% Cell 16 — Print final ADC results

println("\nADC results:")
for i in eachindex(TE_vals)
    println("TE = $(TE_vals[i]) ms")
    println("   Inside  ADC = $(inside_adc[i])   intercept = $(inside_intercepts[i])")
    println("   Outside ADC = $(outside_adc[i])   intercept = $(outside_intercepts[i])")
end


# %% Cell 17 — End

println("\nSignal-vs-b/TE and ADC analysis completed successfully.")
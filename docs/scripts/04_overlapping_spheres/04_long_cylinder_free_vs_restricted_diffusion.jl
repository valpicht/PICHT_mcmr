# ============================================================
# 04_long_cylinder_free_vs_restricted_diffusion.jl
#
# Goal:
# Validate the directional diffusion behavior of an overlapping-sphere
# substrate using a long cylinder-like SWC geometry.
#
# The cylinder is aligned along the x-axis.
#
# Expected behavior:
#   - diffusion along x should be less restricted;
#   - diffusion along y and z should be more restricted;
#   - therefore, diffusion weighting along x should produce stronger
#     signal attenuation and a higher apparent ADC than along y/z.
#
# This script is inspired by the original directional validation script
# used during the project.
#
# Run from the repository root with:
#
#   julia +1.11 --project=. docs/scripts/04_overlapping_spheres/04_long_cylinder_free_vs_restricted_diffusion.jl
# ============================================================

using MCMRSimulator
using MRIBuilder
using CairoMakie
using DelimitedFiles
using Random
using Statistics

include(joinpath(@__DIR__, "utils_overlapping.jl"))

# ============================================================
# Parameters
# ============================================================

const SWC_NAME = "cylinder_soma_neurons_list_new"

const RANDOM_SEED = 1234

# Keep moderate for a fast validation.
const NSPINS = 1_000

# A small b-value set is enough for a directional validation and ADC fit.
const BVALS = [0.0, 0.5, 1.0]

const TE = 40
const TR = 150
const DIFFUSIVITY = 1.0
const BBOX_MARGIN = 2.0

# Manual timestep used in the original validation approach.
#
# Scientific note:
# This is an approximation. It avoids extremely small automatic timesteps
# for small spheres and makes the validation run feasible.
const MANUAL_TIMESTEP = 0.01  # ms

const DIFFUSION_DIRECTIONS = [
    ("x", [1.0, 0.0, 0.0]),
    ("y", [0.0, 1.0, 0.0]),
    ("z", [0.0, 0.0, 1.0]),
]

# ============================================================
# Output paths
# ============================================================

const OUTPUT_DIR = joinpath("docs", "data", "processed")
const FIGURE_DIR = joinpath("docs", "figures")

const SIGNAL_CSV = joinpath(
    OUTPUT_DIR,
    "long_cylinder_directional_signal.csv",
)

const SUMMARY_CSV = joinpath(
    OUTPUT_DIR,
    "long_cylinder_directional_ADC_summary.csv",
)

const SIGNAL_FIG = joinpath(
    FIGURE_DIR,
    "long_cylinder_directional_signal_vs_b.png",
)

const LOG_SIGNAL_FIG = joinpath(
    FIGURE_DIR,
    "long_cylinder_directional_ln_signal_fit.png",
)

const ADC_FIG = joinpath(
    FIGURE_DIR,
    "long_cylinder_directional_ADC_summary.png",
)

# ============================================================
# Sequence and simulation helpers
# ============================================================

"""
    get_directional_sequences(bvals, TE, direction; TR = 150)

Create a set of DWI sequences for one diffusion direction.

A reference sequence is created at the maximum b-value, and lower
b-values are obtained by rescaling the diffusion gradients with
`MRIBuilder.adjust`.

This follows the logic used in the original validation script.
"""
function get_directional_sequences(bvals, TE, direction; TR = 150)
    if maximum(bvals) == 0
        error("maximum(BVALS) must be larger than 0.")
    end

    ref_sequence = DWI(
        bval = maximum(bvals),
        TE = TE,
        TR = TR,
        gradient = (
            orientation = direction,
        ),
        scanner = MRIBuilder.Siemens_Prisma,
    )

    return MRIBuilder.adjust(
        ref_sequence,
        diffusion = (
            scale = sqrt.(bvals ./ maximum(bvals)),
        ),
        merge = false,
    )
end

"""
    simulate_inside_signal_vs_b(spins, geometry, bvals, TE, direction)

Run one combined simulation for all b-values in a given diffusion
direction.

Only inside spins are passed to this function, because the aim is to
validate directional restriction inside the cylinder-like substrate.
"""
function simulate_inside_signal_vs_b(spins, geometry, bvals, TE, direction)
    println("Running one combined simulation for direction = ", direction)

    sequences = get_directional_sequences(
        bvals,
        TE,
        direction;
        TR = TR,
    )

    simulation = Simulation(
        sequences,
        diffusivity = DIFFUSIVITY,
        geometry = geometry,
        timestep = MANUAL_TIMESTEP,
    )

    sigs = readout(
        spins,
        simulation,
        skip_TR = 0,
    )

    signals = [abs(sig.orient.transverse) for sig in sigs]

    return signals
end

"""
    log_normalized_signal(signals)

Compute ln(S/S0), safely handling missing or non-positive values.
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

"""
    fit_adc_from_ln_signal(bvals, ln_signal)

Fit the linear model:

    ln(S/S0) = intercept - ADC * b

using a simple least-squares line fit.
"""
function fit_adc_from_ln_signal(bvals, ln_signal)
    valid = .!isnan.(bvals) .& .!isnan.(ln_signal)

    b = bvals[valid]
    y = ln_signal[valid]

    if length(b) < 2
        return NaN, NaN
    end

    b_mean = mean(b)
    y_mean = mean(y)

    slope = sum((b .- b_mean) .* (y .- y_mean)) / sum((b .- b_mean).^2)
    intercept = y_mean - slope * b_mean

    adc = -slope

    return adc, intercept
end

"""
    fitted_ln_signal(bvals, adc, intercept)

Return fitted values for:

    ln(S/S0) = intercept - ADC * b
"""
function fitted_ln_signal(bvals, adc, intercept)
    return intercept .- adc .* bvals
end

# ============================================================
# Main script
# ============================================================

function main()
    println("============================================================")
    println("04 — Long cylinder: free vs restricted diffusion")
    println("============================================================")

    mkpath(OUTPUT_DIR)
    mkpath(FIGURE_DIR)

    Random.seed!(RANDOM_SEED)

    println("\nSelected SWC substrate:")
    println(SWC_NAME)

    println("\nValidation parameters")
    println("---------------------")
    println("random seed = ", RANDOM_SEED)
    println("nspins = ", NSPINS)
    println("b-values = ", BVALS)
    println("TE = ", TE, " ms")
    println("TR = ", TR, " ms")
    println("diffusivity = ", DIFFUSIVITY)
    println("manual timestep = ", MANUAL_TIMESTEP, " ms")

    geometry, positions, radii = load_overlapping_spheres(SWC_NAME)

    println("\nSubstrate extent:")
    print_sphere_ranges(positions, radii)

    bbox = make_bbox_from_spheres(
        positions,
        radii;
        margin = BBOX_MARGIN,
    )

    println("\nBounding box:")
    println(bbox)

    all_spins = Snapshot(NSPINS, bbox)

    inside_mask = isinside(geometry, all_spins) .> 0
    inside_spins = all_spins[inside_mask]
    outside_spins = all_spins[.!inside_mask]

    println("\nSpin counts")
    println("-----------")
    println("Total spins: ", length(all_spins))
    println("Inside spins: ", length(inside_spins))
    println("Outside spins: ", length(outside_spins))

    if length(inside_spins) == 0
        error("No inside spins were generated. Increase NSPINS or check the geometry.")
    end

    signal_rows = Any[]
    push!(
        signal_rows,
        [
            "TE",
            "direction_name",
            "gradient_x",
            "gradient_y",
            "gradient_z",
            "b",
            "manual_timestep_ms",
            "n_inside_spins",
            "inside_signal",
            "inside_normalized_signal",
            "inside_ln_S_over_S0",
            "inside_fitted_ln_S_over_S0",
        ],
    )

    summary_rows = Any[]
    push!(
        summary_rows,
        [
            "TE",
            "direction_name",
            "gradient_x",
            "gradient_y",
            "gradient_z",
            "manual_timestep_ms",
            "n_inside_spins",
            "inside_ADC",
            "inside_intercept",
        ],
    )

    normalized_by_direction = Dict{String, Vector{Float64}}()
    ln_by_direction = Dict{String, Vector{Float64}}()
    fit_by_direction = Dict{String, Vector{Float64}}()
    adc_by_direction = Dict{String, Float64}()
    intercept_by_direction = Dict{String, Float64}()

    for (direction_name, direction) in DIFFUSION_DIRECTIONS
        println("\n------------------------------------------------------------")
        println("Running direction = ", direction_name)
        println("Gradient orientation = ", direction)
        println("------------------------------------------------------------")

        inside_signals = @time simulate_inside_signal_vs_b(
            inside_spins,
            geometry,
            BVALS,
            TE,
            direction,
        )

        inside_norm = normalize_by_first_value(inside_signals)
        inside_ln = log_normalized_signal(inside_signals)

        inside_adc, inside_intercept = fit_adc_from_ln_signal(
            BVALS,
            inside_ln,
        )

        inside_fit_ln = fitted_ln_signal(
            BVALS,
            inside_adc,
            inside_intercept,
        )

        normalized_by_direction[direction_name] = inside_norm
        ln_by_direction[direction_name] = inside_ln
        fit_by_direction[direction_name] = inside_fit_ln
        adc_by_direction[direction_name] = inside_adc
        intercept_by_direction[direction_name] = inside_intercept

        println("Inside ADC for direction ", direction_name, " = ", inside_adc)
        println("Inside intercept for direction ", direction_name, " = ", inside_intercept)

        for i in eachindex(BVALS)
            push!(
                signal_rows,
                [
                    TE,
                    direction_name,
                    direction[1],
                    direction[2],
                    direction[3],
                    BVALS[i],
                    MANUAL_TIMESTEP,
                    length(inside_spins),
                    inside_signals[i],
                    inside_norm[i],
                    inside_ln[i],
                    inside_fit_ln[i],
                ],
            )
        end

        push!(
            summary_rows,
            [
                TE,
                direction_name,
                direction[1],
                direction[2],
                direction[3],
                MANUAL_TIMESTEP,
                length(inside_spins),
                inside_adc,
                inside_intercept,
            ],
        )

        # Save progressively after each direction.
        open(SIGNAL_CSV, "w") do io
            writedlm(io, signal_rows, ',')
        end

        open(SUMMARY_CSV, "w") do io
            writedlm(io, summary_rows, ',')
        end

        println("Saved intermediate results after direction = ", direction_name)
    end

    println("\nFinal data saved:")
    println(" - ", SIGNAL_CSV)
    println(" - ", SUMMARY_CSV)

    # ========================================================
    # Plot 1 — Normalized signal S(b)/S(0)
    # ========================================================

    println("\nPlotting normalized directional signal...")

    fig_signal = Figure(size = (900, 600))
    ax_signal = Axis(
        fig_signal[1, 1],
        title = "Long overlapping-sphere cylinder: directional diffusion",
        xlabel = "b",
        ylabel = "S(b)/S(0)",
    )

    for (direction_name, _) in DIFFUSION_DIRECTIONS
        lines!(
            ax_signal,
            BVALS,
            normalized_by_direction[direction_name],
            label = direction_name,
        )

        scatter!(
            ax_signal,
            BVALS,
            normalized_by_direction[direction_name],
        )
    end

    axislegend(ax_signal, position = :rb)

    save(SIGNAL_FIG, fig_signal)
    display(fig_signal)

    println("Signal figure saved:")
    println(" - ", SIGNAL_FIG)

    # ========================================================
    # Plot 2 — Log-normalized signal ln(S/S0) with ADC fits
    # ========================================================

    println("\nPlotting log-normalized signal with ADC fits...")

    fig_log = Figure(size = (900, 600))
    ax_log = Axis(
        fig_log[1, 1],
        title = "Directional ADC fit in long overlapping-sphere cylinder",
        xlabel = "b",
        ylabel = "ln(S(b)/S(0))",
    )

    for (direction_name, _) in DIFFUSION_DIRECTIONS
        scatter!(
            ax_log,
            BVALS,
            ln_by_direction[direction_name],
            label = "data " * direction_name,
        )

        lines!(
            ax_log,
            BVALS,
            fit_by_direction[direction_name],
            label = "fit " * direction_name,
        )
    end

    axislegend(ax_log, position = :lb)

    save(LOG_SIGNAL_FIG, fig_log)
    display(fig_log)

    println("Log-signal fit figure saved:")
    println(" - ", LOG_SIGNAL_FIG)

    # ========================================================
    # Plot 3 — ADC summary
    # ========================================================

    println("\nPlotting ADC summary...")

    direction_names = [direction_name for (direction_name, _) in DIFFUSION_DIRECTIONS]
    adc_values = [adc_by_direction[name] for name in direction_names]

    fig_adc = Figure(size = (700, 500))
    ax_adc = Axis(
        fig_adc[1, 1],
        title = "Directional ADC in long overlapping-sphere cylinder",
        xlabel = "Diffusion direction",
        ylabel = "ADC",
        xticks = (1:length(direction_names), direction_names),
    )

    barplot!(ax_adc, 1:length(direction_names), adc_values)

    save(ADC_FIG, fig_adc)
    display(fig_adc)

    println("ADC figure saved:")
    println(" - ", ADC_FIG)

    # ========================================================
    # Final interpretation
    # ========================================================

    println("\nInterpretation")
    println("--------------")
    println("The long cylinder is aligned with the x-axis.")
    println("Diffusion along x is expected to be less restricted than diffusion")
    println("along y and z. Therefore, the x-direction should show stronger")
    println("diffusion-weighted signal attenuation and a higher ADC.")
    println("")
    println("The ln(S/S0) plot makes the ADC interpretation explicit:")
    println("a steeper negative slope corresponds to a higher ADC and therefore")
    println("less restricted diffusion along that direction.")
    println("")
    println("Expected qualitative result:")
    println("    ADC_x > ADC_y and ADC_x > ADC_z")
    println("")
    println("Observed ADC values:")
    for name in direction_names
        println("    ADC_", name, " = ", adc_by_direction[name])
    end

    println("\nLong-cylinder directional validation completed successfully.")
end

main()
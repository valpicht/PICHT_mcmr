# ============================================================
# 01_create_overlapping_spheres.jl
#
# Goal:
# Create an overlapping-sphere geometry for MCMR from an SWC substrate.
#
# This script does not run an MRI simulation yet. It only checks that:
#   1. the SWC file can be read;
#   2. sphere positions and radii can be extracted;
#   3. an MCMR Spheres geometry with overlapping=true can be created;
#   4. a bounding box and random spins can be generated.
#
# Run from the repository root with:
#
#     julia +1.11 --project=. docs/scripts/04_overlapping_spheres/01_create_overlapping_spheres.jl
# ============================================================

using MCMRSimulator
using MRIBuilder

include(joinpath(@__DIR__, "utils_overlapping.jl"))

# ============================================================
# Parameters
# ============================================================

const SWC_NAME = "one_soma_list_new"
const NSPINS = 10_000
const BBOX_MARGIN = 2.0

# ============================================================
# Main script
# ============================================================

function main()
    println("============================================================")
    println("01 — Create overlapping spheres from SWC")
    println("============================================================")

    println("\nSelected SWC substrate:")
    println(SWC_NAME)

    geometry, positions, radii = load_overlapping_spheres(SWC_NAME)

    println("\nSphere coordinate ranges:")
    print_sphere_ranges(positions, radii)

    println("\nCreating random spins around the overlapping-sphere substrate...")
    spins = create_spins_for_spheres(
        positions,
        radii;
        nspins = NSPINS,
        margin = BBOX_MARGIN,
    )

    println("\nCounting inside/outside spins...")
    n_inside, n_outside = print_inside_outside_counts(
        geometry,
        spins,
    )

    println("\nSummary")
    println("-------")
    println("Number of spheres: ", length(radii))
    println("Number of generated spins: ", length(spins))
    println("Inside spins: ", n_inside)
    println("Outside spins: ", n_outside)

    println("\nOverlapping-sphere loading test completed successfully.")
end

main()
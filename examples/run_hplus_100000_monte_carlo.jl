# Weighted Monte Carlo transport of 100,000 precipitating solar-wind protons.
#
# This file is deliberately verbose because it is intended to be inspected and
# modified as a scientific example, rather than used as a hidden production
# driver. Run it from the repository root so that Julia finds Project.toml:
#
#   julia --project=. -t auto examples/run_hplus_100000_monte_carlo.jl
#
# The physical source imposed at the 600 km injection boundary is:
#   species                 H+
#   bulk velocity           [-400, 0, 0] km/s in MSO
#   physical temperature    10 eV
#   physical number density 5 cm^-3 = 5e6 m^-3
#
# The velocity sampler deliberately uses 50 eV, five times the physical
# temperature. This importance-sampling distribution places more numerical
# particles in the velocity tails. MonteCarloWeight corrects every trajectory
# by f(v; 10 eV)/fs(v; 50 eV), so the weighted ensemble still represents the
# requested 10 eV solar-wind proton distribution.
#
# The compact output is a three-dimensional histogram with Julia dimensions
#
#   altitude surface × energy bin × charge state
#
# where charge index 1 is neutral H ENA and charge index 2 is H+. An H ENA
# contribution therefore records charge exchange of an injected proton. Each
# crossing adds the particle density weight once. Upward and downward
# crossings are currently combined, and a recrossing is counted again. The
# result is a crossing-weight diagnostic, not a formally normalized local
# phase-space density.

using MarsASPEN
using MAT
using Printf

# Resolve every path relative to this script. This keeps the example portable
# when the repository is cloned into a different parent directory.
repo = normpath(joinpath(@__DIR__, ".."))
output = length(ARGS) >= 1 ? ARGS[1] :
    joinpath(repo, "examples", "output", "hplus_100000_monte_carlo.mat")

# Load the solar-minimum MGITM neutral atmosphere, extrapolated down to 80 km,
# the MAMPS hot-O component, collision cross sections, energy-loss tables, and
# scattering-angle lookup tables bundled with MarsASPEN.
model = load_model(repo; solar="solar_min", ls=0)

# MonteCarloConfig controls transport and the physical velocity distribution.
# initial_charge_state=1 selects H+. The bulk direction is toward negative MSO
# X in the current one-boundary setup. The seed makes the threaded calculation
# reproducible because each particle receives a deterministic random stream.
config = MonteCarloConfig(
    n_particles=100_000,
    initial_altitude_km=600.0,
    initial_speed_m_s=400_000.0,
    initial_charge_state=1,
    initial_temperature_ev=10.0,
    seed=22,
)

# Source density and importance sampling belong in MonteCarloWeight rather
# than MonteCarloConfig. The SI density 5e6 m^-3 is exactly 5 particles cm^-3.
weighting = MonteCarloWeight(
    sampling_temperature_factor=5.0,
    source_number_density_m3=5.0e6,
)

# Altitude is sampled at the centers of 1 km layers. The final surface is
# 599.5 km, immediately below the 600 km source boundary. Energy uses 100
# logarithmic bins between 1 and 10,000 eV. The plotting example displays only
# 10 to 3,000 eV, but retaining the wider grid helps diagnose outliers.
altitude_edges_km = collect(80.0:1.0:600.0)
altitude_centers_km = 0.5 .* (
    altitude_edges_km[1:end-1] .+ altitude_edges_km[2:end]
)
energy_edges_ev = 10.0 .^ range(0.0, 4.0, length=101)

# Only the compact crossing histogram is retained. Saving every event from
# 100,000 trajectories would create a much larger MAT file and is unnecessary
# for the altitude-energy figure.
elapsed = @elapsed result = run_density_crossing_ensemble(
    model, config;
    weighting=weighting,
    altitude_surfaces_km=altitude_centers_km,
    energy_edges_ev=energy_edges_ev,
)

# MAT v7.3 is used by MAT.jl for this multidimensional result. The companion
# Python reader supports both ordinary MAT files and HDF5-backed MAT v7.3.
mkpath(dirname(abspath(output)))
matwrite(output, Dict(
    "format_version" => "marsaspen_hplus_weighted_example_v1",
    "n_particles" => config.n_particles,
    "seed" => config.seed,
    "initial_species" => "Hplus",
    "initial_altitude_km" => config.initial_altitude_km,
    "initial_bulk_velocity_m_s" =>
        [-config.initial_speed_m_s, 0.0, 0.0],
    "physical_temperature_ev" => config.initial_temperature_ev,
    "sampling_temperature_ev" =>
        config.initial_temperature_ev *
        weighting.sampling_temperature_factor,
    "sampling_temperature_factor" =>
        weighting.sampling_temperature_factor,
    "source_number_density_m3" =>
        weighting.source_number_density_m3,
    "total_importance_weight" => result.total_importance_weight,
    "altitude_edges_km" => altitude_edges_km,
    "energy_edges_ev" => energy_edges_ev,
    "charge_state_names" => ["H_ENA", "Hplus"],
    "density_weight_sum_m3" => result.density_weight_sum_m3,
))

# These diagnostics make performance and importance-weight convergence visible
# in terminal logs. For a sufficiently large ensemble, total_importance_weight
# should be close to n_particles.
@printf("n_particles=%d\n", config.n_particles)
@printf("elapsed_s=%.6f\n", elapsed)
@printf("particles_per_s=%.3f\n", config.n_particles / elapsed)
@printf("total_importance_weight=%.9g\n", result.total_importance_weight)
@printf("output=%s\n", abspath(output))

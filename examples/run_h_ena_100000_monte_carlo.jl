# Weighted Monte Carlo example for 100,000 precipitating H ENA particles.
#
# Run from the repository root:
#   julia --project=. -t auto examples/run_h_ena_100000_monte_carlo.jl
#
# Physical source:
#   initial altitude       600 km
#   bulk velocity          [-400, 0, 0] km/s in MSO
#   physical temperature   10 eV
#   source number density  1 cm^-3 = 1e6 m^-3
#
# Importance sampling:
#   Velocities are sampled from a 50 eV drifting Maxwellian, which is five
#   times broader than the physical distribution. Each trajectory receives
#   weight f(v; 10 eV) / fs(v; 50 eV), normalized so the density weights sum
#   exactly to 1e6 m^-3.
#
# Output:
#   A compact MAT file containing the sum of particle density weights in each
#   altitude and energy bin. No path-length factor and no division by energy
#   width are applied. Full trajectories are intentionally not saved.
#
# Array layout is altitude surface × energy bin × charge state. Charge index 1
# is H ENA and index 2 is H+. Upward and downward crossings are combined, and
# a returning particle can cross a surface more than once.

using MarsASPEN
using MAT
using Printf

# Resolve paths relative to the cloned repository.
repo = normpath(joinpath(@__DIR__, ".."))
output = length(ARGS) >= 1 ? ARGS[1] :
    joinpath(repo, "examples", "output", "h_ena_100000_monte_carlo.mat")

# Assemble the solar-minimum atmosphere, hot O, collision cross sections,
# reaction energy losses, and scattering-angle data.
model = load_model(repo; solar="solar_min", ls=0)

# Transport settings contain no source-density or macro-particle weighting.
# State zero selects neutral H. The 10 eV temperature broadens velocities about
# the 400 km/s bulk flow. The seed gives deterministic particle random streams.
config = MonteCarloConfig(
    n_particles=100_000,
    initial_altitude_km=600.0,
    initial_speed_m_s=400_000.0,
    initial_charge_state=0,
    initial_temperature_ev=10.0,
    seed=21,
)

# The 50 eV sampling distribution improves tail statistics. The f/fs weight
# restores the physical 10 eV source, and 1 cm^-3 is supplied as 1e6 m^-3.
weighting = MonteCarloWeight(
    sampling_temperature_factor=5.0,
    source_number_density_m3=1.0e6,
)

# One-kilometer altitude bins and 100 logarithmic energy bins cover 1 to
# 10,000 eV. A logarithmic grid cannot begin at exactly zero. Because transport
# stops below 10 eV, the lowest part of this requested range is normally empty.
altitude_edges_km = collect(80.0:1.0:600.0)
altitude_centers_km = 0.5 .* (
    altitude_edges_km[1:end-1] .+ altitude_edges_km[2:end]
)
energy_edges_ev = 10.0 .^ range(0.0, 4.0, length=101)

# Compact mode retains only the histogram, not every particle event.
elapsed = @elapsed result = run_density_crossing_ensemble(
    model, config;
    weighting=weighting,
    altitude_surfaces_km=altitude_centers_km,
    energy_edges_ev=energy_edges_ev,
)

# Store source metadata with the histogram so Python does not infer parameters
# from the filename. The reader also corrects reversed MAT v7.3 dimensions.
mkpath(dirname(abspath(output)))
matwrite(output, Dict(
    "format_version" => "marsaspen_h_ena_weighted_example_v1",
    "n_particles" => config.n_particles,
    "seed" => config.seed,
    "initial_species" => "H_ENA",
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

# Print performance and importance-weight convergence diagnostics.
@printf("n_particles=%d\n", config.n_particles)
@printf("elapsed_s=%.6f\n", elapsed)
@printf("particles_per_s=%.3f\n", config.n_particles / elapsed)
@printf("total_importance_weight=%.9g\n", result.total_importance_weight)
@printf("output=%s\n", abspath(output))

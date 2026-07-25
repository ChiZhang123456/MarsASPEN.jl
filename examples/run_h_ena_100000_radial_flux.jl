# Local radial-flux profile for 100,000 precipitating H ENA particles.
#
# The physical source is placed at the 600 km spherical boundary:
#
#   initial species         neutral H ENA
#   bulk velocity           400 km/s radially inward
#   physical temperature    10 eV
#   physical number density 5 cm^-3 = 5e6 m^-3
#
# Every numerical particle has a normalized density weight Wn in m^-3. At
# every altitude-surface crossing, MarsASPEN evaluates the current local radial
# velocity Vr after previous scattering, energy loss, and charge-state changes.
# It then accumulates Wn*abs(Vr), in m^-2 s^-1, into a downward or upward bin.
#
# Although the injected population is neutral, stripping can create H+ during
# transport. The saved arrays therefore contain both H ENA and H+ profiles.
#
# Run from the repository root:
#
#   julia --project=. -t auto examples/run_h_ena_100000_radial_flux.jl

using MarsASPEN
using MAT
using Printf

repo = normpath(joinpath(@__DIR__, ".."))
output = length(ARGS) >= 1 ? ARGS[1] :
    joinpath(repo, "examples", "output", "h_ena_100000_radial_flux.mat")

model = load_model(repo; solar="solar_min", ls=0)
config = MonteCarloConfig(
    n_particles=100_000,
    initial_altitude_km=600.0,
    initial_speed_m_s=400_000.0,
    initial_charge_state=0,
    initial_temperature_ev=10.0,
    seed=23,
)

# A broader sampling distribution improves coverage of the thermal tails.
# Importance weighting restores the requested 10 eV physical distribution,
# and the normalized particle density weights sum to 5e6 m^-3.
weighting = MonteCarloWeight(
    sampling_temperature_factor=5.0,
    source_number_density_m3=5.0e6,
)

altitude_surfaces_km = collect(80.5:1.0:599.5)
elapsed = @elapsed result = run_radial_flux_ensemble(
    model, config;
    weighting=weighting,
    altitude_surfaces_km=altitude_surfaces_km,
)

mkpath(dirname(abspath(output)))
matwrite(output, Dict(
    "format_version" => "marsaspen_local_radial_flux_v1",
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
    "source_number_density_m3" => weighting.source_number_density_m3,
    "nominal_bulk_flux_m2_s" =>
        weighting.source_number_density_m3 * config.initial_speed_m_s,
    "total_importance_weight" => result.total_importance_weight,
    "altitude_surfaces_km" => result.altitude_surfaces_km,
    "charge_state_names" => ["H_ENA", "Hplus"],
    "direction_names" => ["downward", "upward"],
    "radial_flux_m2_s" => result.radial_flux_m2_s,
    "signed_outward_flux_m2_s" => result.signed_outward_flux_m2_s,
    "net_downward_flux_m2_s" => result.net_downward_flux_m2_s,
))

@printf("n_particles=%d\n", config.n_particles)
@printf("elapsed_s=%.6f\n", elapsed)
@printf("particles_per_s=%.3f\n", config.n_particles / elapsed)
@printf("nominal_bulk_flux_m2_s=%.9g\n",
        weighting.source_number_density_m3 * config.initial_speed_m_s)
@printf("top_HENA_downward_flux_m2_s=%.9g\n",
        result.radial_flux_m2_s[end, 1, 1])
@printf("output=%s\n", abspath(output))

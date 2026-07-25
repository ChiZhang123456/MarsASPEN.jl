# Geographic diagnostic maps at the 120 km spherical surface.
#
# Usage:
#
#   julia --project=. -t auto examples/run_surface_diagnostics_120km.jl h_ena
#   julia --project=. -t auto examples/run_surface_diagnostics_120km.jl hplus
#
# This first implementation deliberately retains the original MarsASPEN point
# source at lon=0 deg, lat=0 deg and 600 km altitude. The output is therefore a
# local precipitation footprint, not a globally normalized dayside map.
#
# Ionization and Ly-alpha rates use the local total projectile speed |V|.
# The collision-energy field includes expected kinetic-energy transfer by all
# modeled reaction channels. The terminal field additionally deposits the
# remaining energy of particles stopped below 10 eV within 119.5 to 120.5 km.

using MarsASPEN
using MAT
using Printf

repo = normpath(joinpath(@__DIR__, ".."))
source_name = length(ARGS) >= 1 ? lowercase(ARGS[1]) : "h_ena"
source_name in ("h_ena", "hplus") ||
    throw(ArgumentError("source must be h_ena or hplus"))
initial_charge = source_name == "hplus" ? 1 : 0
output = length(ARGS) >= 2 ? ARGS[2] : joinpath(
    repo, "examples", "output",
    "$(source_name)_100000_surface_diagnostics_120km.mat",
)

model = load_model(repo; solar="solar_min", ls=0)
config = MonteCarloConfig(
    n_particles=100_000,
    initial_altitude_km=600.0,
    initial_speed_m_s=400_000.0,
    initial_charge_state=initial_charge,
    initial_temperature_ev=10.0,
    seed=initial_charge == 1 ? 52 : 51,
    min_energy_ev=10.0,
    include_hot_o=true,
)
weighting = MonteCarloWeight(
    sampling_temperature_factor=5.0,
    source_number_density_m3=5.0e6,
)

elapsed = @elapsed result = run_surface_diagnostic_map_ensemble(
    model, config;
    weighting=weighting,
    altitude_km=120.0,
    longitude_edges_deg=collect(-30.0:1.0:30.0),
    latitude_edges_deg=collect(-30.0:1.0:30.0),
    thermalization_shell_half_width_km=0.5,
)

mkpath(dirname(abspath(output)))
matwrite(output, Dict(
    "format_version" => "marsaspen_point_source_surface_map_v1",
    "source_geometry" => result.source_geometry,
    "initial_species" => initial_charge == 1 ? "Hplus" : "H_ENA",
    "n_particles" => config.n_particles,
    "seed" => config.seed,
    "initial_altitude_km" => config.initial_altitude_km,
    "initial_speed_m_s" => config.initial_speed_m_s,
    "physical_temperature_ev" => config.initial_temperature_ev,
    "source_number_density_m3" => weighting.source_number_density_m3,
    "altitude_km" => result.altitude_km,
    "longitude_edges_deg" => result.longitude_edges_deg,
    "latitude_edges_deg" => result.latitude_edges_deg,
    "oxygen_ionization_rate_m3_s1" =>
        result.oxygen_ionization_rate_m3_s1,
    "co2_ionization_rate_m3_s1" =>
        result.co2_ionization_rate_m3_s1,
    "total_ionization_rate_m3_s1" =>
        result.total_ionization_rate_m3_s1,
    "lya_volume_emission_rate_photons_m3_s1" =>
        result.lya_volume_emission_rate_photons_m3_s1,
    "collision_energy_transfer_ev_m3_s1" =>
        result.collision_energy_transfer_ev_m3_s1,
    "thermalized_below_cutoff_ev_m3_s1" =>
        result.thermalized_below_cutoff_ev_m3_s1,
    "total_energy_transfer_ev_m3_s1" =>
        result.total_energy_transfer_ev_m3_s1,
    "total_energy_transfer_w_m3" =>
        result.total_energy_transfer_w_m3,
    "thermalization_shell_edges_km" =>
        result.thermalization_shell_edges_km,
    "elastic_mean_loss_fraction" =>
        result.elastic_mean_loss_fraction,
))

@printf("source=%s\n", source_name)
@printf("elapsed_s=%.6f\n", elapsed)
@printf("max_total_ionization_rate_m3_s1=%.9g\n",
        maximum(result.total_ionization_rate_m3_s1))
@printf("max_lya_VER_photons_m3_s1=%.9g\n",
        maximum(result.lya_volume_emission_rate_photons_m3_s1))
@printf("max_total_energy_transfer_w_m3=%.9g\n",
        maximum(result.total_energy_transfer_w_m3))
@printf("output=%s\n", abspath(output))

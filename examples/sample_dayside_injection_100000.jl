# Sample the global dayside solar-wind injection boundary without transport.
#
# The 100,000 initial H+ macro particles are uniformly distributed in surface
# area over the 600 km dayside hemisphere. Their MSO velocity distribution is
# a drifting 3D Maxwellian with bulk velocity (-400, 0, 0) km/s and kT=10 eV.
# Local radial velocity is evaluated independently at every position.
#
# Usage:
#
#   julia --project=. -t auto examples/sample_dayside_injection_100000.jl

using MarsASPEN
using MAT
using Printf
using Statistics

repo = normpath(joinpath(@__DIR__, ".."))
output = length(ARGS) >= 1 ? ARGS[1] : joinpath(
    repo, "examples", "output", "dayside_hplus_injection_100000.mat",
)

config = MonteCarloConfig(
    n_particles=100_000,
    injection_geometry=:dayside_uniform,
    initial_altitude_km=600.0,
    initial_speed_m_s=400_000.0,
    initial_charge_state=1,
    initial_temperature_ev=10.0,
    seed=61,
)
weighting = MonteCarloWeight(
    sampling_temperature_factor=1.0,
    source_number_density_m3=5.0e6,
)

sample = sample_injection_ensemble(config; weighting=weighting)
velocity = sample.velocity_m_s

mkpath(dirname(abspath(output)))
matwrite(output, Dict(
    "format_version" => "marsaspen_dayside_injection_v1",
    "n_particles" => config.n_particles,
    "seed" => config.seed,
    "injection_geometry" => sample.injection_geometry,
    "initial_species" => "Hplus",
    "initial_altitude_km" => config.initial_altitude_km,
    "bulk_velocity_m_s" => [-config.initial_speed_m_s, 0.0, 0.0],
    "physical_temperature_ev" => config.initial_temperature_ev,
    "source_number_density_m3" => weighting.source_number_density_m3,
    "position_m" => sample.position_m,
    "velocity_m_s" => velocity,
    "longitude_deg" => sample.longitude_deg,
    "latitude_deg" => sample.latitude_deg,
    "solar_zenith_angle_deg" => sample.solar_zenith_angle_deg,
    "radial_velocity_m_s" => sample.radial_velocity_m_s,
    "density_weight_m3" => sample.density_weight_m3,
    "inward_flux_weight_m2_s" => sample.inward_flux_weight_m2_s,
    "importance_weight" => sample.importance_weight,
))

@printf("n_particles=%d\n", config.n_particles)
@printf("mean_velocity_km_s=(%.6f, %.6f, %.6f)\n",
        mean(velocity[:, 1]) / 1000,
        mean(velocity[:, 2]) / 1000,
        mean(velocity[:, 3]) / 1000)
@printf("std_velocity_km_s=(%.6f, %.6f, %.6f)\n",
        std(velocity[:, 1]) / 1000,
        std(velocity[:, 2]) / 1000,
        std(velocity[:, 3]) / 1000)
@printf("SZA_range_deg=(%.6f, %.6f)\n",
        minimum(sample.solar_zenith_angle_deg),
        maximum(sample.solar_zenith_angle_deg))
@printf("inward_fraction=%.9f\n",
        mean(sample.radial_velocity_m_s .< 0))
@printf("sum_density_weight_m3=%.9g\n", sum(sample.density_weight_m3))
@printf("output=%s\n", abspath(output))

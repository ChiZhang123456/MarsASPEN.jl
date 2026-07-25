# H Ly-alpha volume-emission profiles from 100,000 weighted macro particles.
#
# Usage:
#
#   julia --project=. -t auto examples/run_lya_volume_emission.jl h_ena
#   julia --project=. -t auto examples/run_lya_volume_emission.jl hplus
#
# Both source cases use 600 km, 400 km/s inward bulk speed, 10 eV physical
# temperature, 5 cm^-3 source density, solar minimum, and Ls=0. The effective
# Ly-alpha production cross sections are evaluated at each crossing energy for
# H ENA and H+ impacts on CO2, O, and N2.
#
# The saved volume emission rate has units photons m^-3 s^-1. Multiplication by
# h*c/lambda at 121.567 nm gives the radiative source in W m^-3. This quantity
# is radiative energy, not local atmospheric heating.

using MarsASPEN
using MAT
using Printf

repo = normpath(joinpath(@__DIR__, ".."))
source_name = length(ARGS) >= 1 ? lowercase(ARGS[1]) : "h_ena"
source_name in ("h_ena", "hplus") ||
    throw(ArgumentError("source must be h_ena or hplus"))
initial_charge = source_name == "hplus" ? 1 : 0
default_output = joinpath(
    repo, "examples", "output",
    "$(source_name)_100000_lya_volume_emission.mat",
)
output = length(ARGS) >= 2 ? ARGS[2] : default_output

model = load_model(repo; solar="solar_min", ls=0)
config = MonteCarloConfig(
    n_particles=100_000,
    initial_altitude_km=600.0,
    initial_speed_m_s=400_000.0,
    initial_charge_state=initial_charge,
    initial_temperature_ev=10.0,
    seed=initial_charge == 1 ? 42 : 41,
    include_hot_o=true,
)
weighting = MonteCarloWeight(
    sampling_temperature_factor=5.0,
    source_number_density_m3=5.0e6,
)
altitude_surfaces_km = collect(80.5:1.0:599.5)

elapsed = @elapsed result = run_lya_volume_emission_ensemble(
    model, config;
    weighting=weighting,
    altitude_surfaces_km=altitude_surfaces_km,
    density_lon_deg=0.0,
    density_lat_deg=0.0,
)

mkpath(dirname(abspath(output)))
matwrite(output, Dict(
    "format_version" => "marsaspen_lya_volume_emission_v1",
    "n_particles" => config.n_particles,
    "seed" => config.seed,
    "initial_species" => initial_charge == 1 ? "Hplus" : "H_ENA",
    "initial_altitude_km" => config.initial_altitude_km,
    "initial_speed_m_s" => config.initial_speed_m_s,
    "physical_temperature_ev" => config.initial_temperature_ev,
    "source_number_density_m3" => weighting.source_number_density_m3,
    "altitude_surfaces_km" => result.altitude_surfaces_km,
    "charge_state_names" => ["H_ENA", "Hplus"],
    "target_names" => ["CO2", "O", "N2"],
    "direction_names" => ["downward", "upward"],
    "target_density_m3" => result.target_density_m3,
    "lya_rate_coefficient_s1" => result.lya_rate_coefficient_s1,
    "volume_emission_rate_photons_m3_s1" =>
        result.volume_emission_rate_photons_m3_s1,
    "volume_emission_rate_by_charge_photons_m3_s1" =>
        result.volume_emission_rate_by_charge_photons_m3_s1,
    "volume_emission_rate_by_target_photons_m3_s1" =>
        result.volume_emission_rate_by_target_photons_m3_s1,
    "total_volume_emission_rate_photons_m3_s1" =>
        result.total_volume_emission_rate_photons_m3_s1,
    "photon_energy_j" => result.photon_energy_j,
    "radiative_energy_rate_by_charge_w_m3" =>
        result.radiative_energy_rate_by_charge_w_m3,
    "total_radiative_energy_rate_w_m3" =>
        result.total_radiative_energy_rate_w_m3,
))

peak_ver, peak_index =
    findmax(result.total_volume_emission_rate_photons_m3_s1)
@printf("source=%s\n", source_name)
@printf("elapsed_s=%.6f\n", elapsed)
@printf("photon_energy_j=%.9g\n", result.photon_energy_j)
@printf("peak_altitude_km=%.1f\n", result.altitude_surfaces_km[peak_index])
@printf("peak_VER_photons_m3_s1=%.9g\n", peak_ver)
@printf("output=%s\n", abspath(output))

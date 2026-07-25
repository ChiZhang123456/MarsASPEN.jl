# Oxygen ionization-rate profile from a weighted 100,000-particle simulation.
#
# Usage:
#
#   julia --project=. -t auto examples/run_oxygen_ionization_rate.jl h_ena
#   julia --project=. -t auto examples/run_oxygen_ionization_rate.jl hplus
#
# Both cases use a 600 km source, 400 km/s inward bulk speed, 10 eV physical
# temperature, and 5 cm^-3 source density. The only difference is the initial
# projectile state. Collisions can subsequently convert H ENA to H+ or H+ to
# H ENA, so each output retains both projectile contributions.
#
# At every spherical altitude crossing, the estimator accumulates
#
#   Wn_i * abs(Vr_i) * sigma_ion(E_i)             [s^-1]
#
# where E_i is calculated from the local total speed after previous collisions.
# Multiplication by the local O density from MGITM cold O plus MAMPS hot O gives
#
#   q_O                                             [m^-3 s^-1].
#
# This is a crossing estimator over all particles. Selecting only realized
# ionization events and multiplying those events by sigma again would count
# the collision probability twice.

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
    "$(source_name)_100000_oxygen_ionization_rate.mat",
)
output = length(ARGS) >= 2 ? ARGS[2] : default_output

model = load_model(repo; solar="solar_min", ls=0)
config = MonteCarloConfig(
    n_particles=100_000,
    initial_altitude_km=600.0,
    initial_speed_m_s=400_000.0,
    initial_charge_state=initial_charge,
    initial_temperature_ev=10.0,
    seed=initial_charge == 1 ? 32 : 31,
    include_hot_o=true,
)
weighting = MonteCarloWeight(
    sampling_temperature_factor=5.0,
    source_number_density_m3=5.0e6,
)
altitude_surfaces_km = collect(80.5:1.0:599.5)

elapsed = @elapsed result = run_target_ionization_rate_ensemble(
    model, config;
    weighting=weighting,
    altitude_surfaces_km=altitude_surfaces_km,
    target=:O,
    density_lon_deg=0.0,
    density_lat_deg=0.0,
)

mkpath(dirname(abspath(output)))
matwrite(output, Dict(
    "format_version" => "marsaspen_target_ionization_rate_v1",
    "n_particles" => config.n_particles,
    "seed" => config.seed,
    "initial_species" => initial_charge == 1 ? "Hplus" : "H_ENA",
    "initial_altitude_km" => config.initial_altitude_km,
    "initial_speed_m_s" => config.initial_speed_m_s,
    "physical_temperature_ev" => config.initial_temperature_ev,
    "source_number_density_m3" => weighting.source_number_density_m3,
    "target_name" => result.target_name,
    "target_density_m3" => result.target_density_m3,
    "density_lon_deg" => result.density_lon_deg,
    "density_lat_deg" => result.density_lat_deg,
    "altitude_surfaces_km" => result.altitude_surfaces_km,
    "charge_state_names" => ["H_ENA", "Hplus"],
    "direction_names" => ["downward", "upward"],
    "radial_flux_m2_s" => result.radial_flux_m2_s,
    "flux_times_ionization_cross_section_s1" =>
        result.flux_times_ionization_cross_section_s1,
    "ionization_rate_m3_s1" => result.ionization_rate_m3_s1,
    "ionization_rate_by_charge_m3_s1" =>
        result.ionization_rate_by_charge_m3_s1,
    "total_ionization_rate_m3_s1" =>
        result.total_ionization_rate_m3_s1,
))

peak_rate, peak_index = findmax(result.total_ionization_rate_m3_s1)
@printf("source=%s\n", source_name)
@printf("elapsed_s=%.6f\n", elapsed)
@printf("peak_altitude_km=%.1f\n", result.altitude_surfaces_km[peak_index])
@printf("peak_total_O_ionization_rate_m3_s1=%.9g\n", peak_rate)
@printf("output=%s\n", abspath(output))

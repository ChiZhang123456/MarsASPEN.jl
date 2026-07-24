# Solar-wind proton Monte Carlo run with directional flux diagnostics.
#
# Default source at 600 km:
#   H+ number density = 5 cm^-3 = 5e6 m^-3
#   bulk velocity = [-400, 0, 0] km/s
#   drifting-Maxwellian temperature kT = 10 eV
#
# Arguments: particle count, output MAT, altitude spacing in km, energy spacing
# in eV. The output stores upward and downward crossing flux for H ENA and H+.

using MarsASPEN
using MAT
using Printf

repo = normpath(joinpath(@__DIR__, ".."))
n = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 10_000_000
output = length(ARGS) >= 2 ? ARGS[2] :
    joinpath(repo, "output", "solar_wind_flux_$(n)p.mat")
altitude_spacing_km = length(ARGS) >= 3 ? parse(Float64, ARGS[3]) : 0.5
energy_spacing_ev = length(ARGS) >= 4 ? parse(Float64, ARGS[4]) : 2.0

number_density_m3 = 5.0e6
bulk_speed_m_s = 400_000.0
temperature_ev = 10.0
sampling_temperature_factor = 5.0

model = load_model(repo; solar="solar_min", ls=0)
config = MonteCarloConfig(
    n_particles=n,
    initial_altitude_km=600.0,
    initial_speed_m_s=bulk_speed_m_s,
    initial_charge_state=1,
    initial_temperature_ev=temperature_ev,
    sampling_temperature_factor=sampling_temperature_factor,
    initial_number_density_m3=number_density_m3,
    seed=7,
)
altitude_surfaces = collect(100.0:altitude_spacing_km:300.0)
energy_edges = collect(10.0:energy_spacing_ev:2000.0)

# Compile the specialized diagnostic path with small histograms.
run_directional_flux_ensemble(
    model, MonteCarloConfig(
        n_particles=2,
        initial_charge_state=1,
        initial_temperature_ev=temperature_ev,
        sampling_temperature_factor=sampling_temperature_factor,
        initial_number_density_m3=number_density_m3,
    );
    altitude_surfaces_km=collect(100.0:20.0:300.0),
    energy_edges_ev=collect(10.0:100.0:2010.0),
)

elapsed = @elapsed result = run_directional_flux_ensemble(
    model, config;
    altitude_surfaces_km=altitude_surfaces,
    energy_edges_ev=energy_edges,
)

mkpath(dirname(abspath(output)))
matwrite(output, Dict(
    "format_version" => "marsaspen_directional_flux_v1",
    "n_particles" => n,
    "seed" => config.seed,
    "ls_deg" => 0,
    "f107" => 70,
    "initial_species" => "Hplus",
    "initial_altitude_km" => config.initial_altitude_km,
    "initial_bulk_velocity_m_s" => [-bulk_speed_m_s, 0.0, 0.0],
    "initial_number_density_m3" => number_density_m3,
    "initial_temperature_ev" => temperature_ev,
    "sampling_temperature_ev" => temperature_ev * sampling_temperature_factor,
    "sampling_temperature_factor" => sampling_temperature_factor,
    "total_importance_weight" => result.total_importance_weight,
    "nominal_incident_flux_m2_s" => number_density_m3 * bulk_speed_m_s,
    "altitude_surfaces_km" => altitude_surfaces,
    "energy_edges_ev" => energy_edges,
    "charge_state_names" => ["H_ENA", "Hplus"],
    "direction_names" => ["downward", "upward"],
    "flux_m2_s" => result.flux_m2_s,
    "stop_counts" => result.stop_counts,
    "final_energy_mean_ev" => result.final_energy_mean_ev,
    "collision_mean" => result.collision_mean,
    "step_mean" => result.step_mean,
))

@printf("n_particles=%d\n", n)
@printf("threads=%d\n", Threads.nthreads())
@printf("elapsed_s=%.6f\n", elapsed)
@printf("particles_per_s=%.3f\n", n / elapsed)
@printf("nominal_incident_flux_m2_s=%.9g\n", number_density_m3 * bulk_speed_m_s)
@printf("final_energy_mean_ev=%.6f\n", result.final_energy_mean_ev)
@printf("collision_mean=%.6f\n", result.collision_mean)
@printf("step_mean=%.6f\n", result.step_mean)
for code in 1:5
    @printf("stop_%d=%d\n", code, result.stop_counts[code])
end
@printf("output=%s\n", abspath(output))

# Local radial-flux profile for 100,000 precipitating solar-wind protons.
#
# Physical source at the 600 km spherical boundary:
#   initial species         H+
#   bulk velocity           400 km/s radially inward
#   physical temperature    10 eV
#   physical number density 5 cm^-3 = 5e6 m^-3
#
# Each numerical particle carries a density weight
#
#   Wn_i = n_source * (f/fs)_i / sum_j(f/fs)_j        [m^-3].
#
# When a particle crosses an altitude surface, this example evaluates its
# local radial velocity
#
#   Vr = v dot r_hat = (x*vx + y*vy + z*vz) / r       [m s^-1]
#
# and adds Wn_i*abs(Vr) to either the downward or upward flux accumulator.
# Therefore every accumulated flux has units m^-2 s^-1. The calculation uses
# local velocity after all previous energy losses and scattering events. This
# differs from run_directional_flux_ensemble, whose legacy diagnostic uses the
# injection-plane velocity weight.
#
# Run from the repository root:
#
#   julia --project=. -t auto examples/run_hplus_100000_radial_flux.jl

using MarsASPEN
using MAT
using Printf

repo = normpath(joinpath(@__DIR__, ".."))
output = length(ARGS) >= 1 ? ARGS[1] :
    joinpath(repo, "examples", "output", "hplus_100000_radial_flux.mat")

# Use the same solar-minimum, Ls=0 atmosphere as the existing 100,000-particle
# H+ altitude-energy example so the two diagnostics can be compared directly.
model = load_model(repo; solar="solar_min", ls=0)
config = MonteCarloConfig(
    n_particles=100_000,
    initial_altitude_km=600.0,
    initial_speed_m_s=400_000.0,
    initial_charge_state=1,
    initial_temperature_ev=10.0,
    seed=22,
)

# Sample a broader 50 eV Maxwellian and correct each trajectory by f/fs. The
# normalized Wn values sum exactly to the physical 5e6 m^-3 source density.
weighting = MonteCarloWeight(
    sampling_temperature_factor=5.0,
    source_number_density_m3=5.0e6,
)

# One surface is placed at the center of every one-kilometer altitude layer.
# The uppermost surface is 599.5 km, immediately below the source boundary.
altitude_surfaces_km = collect(80.5:1.0:599.5)

elapsed = @elapsed result = run_radial_flux_ensemble(
    model, config;
    weighting=weighting,
    altitude_surfaces_km=altitude_surfaces_km,
)

# Save downward and upward positive magnitudes, plus signed outward and net
# downward profiles. Charge order is H ENA followed by H+.
mkpath(dirname(abspath(output)))
matwrite(output, Dict(
    "format_version" => "marsaspen_local_radial_flux_v1",
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
@printf("top_Hplus_downward_flux_m2_s=%.9g\n",
        result.radial_flux_m2_s[end, 2, 1])
@printf("output=%s\n", abspath(output))

# Full three-dimensional Monte Carlo output for a uniform dayside H+ source.
#
# Usage:
#
#   julia --project=. -t auto examples/run_dayside_3d_10000000.jl
#
# Initial conditions:
#   altitude:            600 km
#   position:            uniform in area over the MSO dayside hemisphere
#   bulk velocity:       (-400, 0, 0) km/s in the global MSO frame
#   physical kT:         10 eV
#   source density:      5 cm^-3 = 5e6 m^-3
#   number of particles: 10,000,000
#
# Horizontal grid edges match the native 5 degree MGITM grid. Altitude edges
# run from 80 to 600 km at 1 km spacing. The example writes four files:
#
#   *_grid.mat       coordinate edges, centers, and spherical cell volumes
#   *_moments.mat    density and total/radial/upward/downward flux
#   *_reactions.mat  state change, ionization, Ly-alpha, and elastic rates
#   *_energy.mat     collision loss and sub-10 eV thermalization
#
# Arrays use dimension order longitude, latitude, altitude, followed by the
# optional charge, target, and reaction dimensions documented in each file.

using MarsASPEN
using Printf

repo = normpath(joinpath(@__DIR__, ".."))
output_prefix = length(ARGS) >= 1 ? ARGS[1] : joinpath(
    repo, "examples", "output", "dayside_hplus_10000000_3d",
)
n_particles = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 10_000_000
minimum_grid_altitude_km =
    length(ARGS) >= 3 ? parse(Float64, ARGS[3]) : 80.0

model = load_model(repo; solar="solar_min", ls=0)
config = MonteCarloConfig(
    n_particles=n_particles,
    injection_geometry=:dayside_uniform,
    initial_altitude_km=600.0,
    initial_speed_m_s=400_000.0,
    initial_charge_state=1,
    initial_temperature_ev=10.0,
    seed=71,
    max_step_m=1000.0,
    min_energy_ev=10.0,
    min_altitude_km=minimum_grid_altitude_km,
    max_altitude_km=1000.0,
    include_hot_o=true,
)
weighting = MonteCarloWeight(
    sampling_temperature_factor=1.0,
    source_number_density_m3=5.0e6,
)

elapsed_transport = @elapsed result = run_spatial_grid_ensemble(
    model, config;
    weighting=weighting,
    altitude_edges_km=collect(minimum_grid_altitude_km:1.0:600.0),
)

metadata = Dict(
    "n_particles" => config.n_particles,
    "seed" => config.seed,
    "initial_species" => "Hplus",
    "injection_geometry" => "uniform dayside spherical surface",
    "coordinate_system" => "MSO",
    "initial_altitude_km" => config.initial_altitude_km,
    "bulk_velocity_m_s" => [-config.initial_speed_m_s, 0.0, 0.0],
    "physical_temperature_ev" => config.initial_temperature_ev,
    "source_number_density_m3" => weighting.source_number_density_m3,
    "dayside_injection_area_m2" => result.dayside_injection_area_m2,
    "stop_counts" => result.stop_counts,
    "total_collisions" => result.total_collisions,
    "total_steps" => result.total_steps,
)
elapsed_write = @elapsed files = write_spatial_grid_mats(
    output_prefix, result.grid; metadata=metadata,
)

@printf("transport_elapsed_s=%.6f\n", elapsed_transport)
@printf("write_elapsed_s=%.6f\n", elapsed_write)
@printf("total_steps=%d\n", result.total_steps)
@printf("total_collisions=%d\n", result.total_collisions)
@printf("grid_file=%s\n", files.grid)
@printf("moments_file=%s\n", files.moments)
@printf("reactions_file=%s\n", files.reactions)
@printf("energy_file=%s\n", files.energy)

# Generate detailed trajectories for an X-R particle-motion animation.
#
# The particles are sampled uniformly over the 600 km dayside hemisphere.
# Their velocities follow the same drifting Maxwellian used by the production
# dayside example. Detailed histories are appropriate here because the
# animation needs the position and charge state of every particle as a
# function of time.
#
# Run from the repository root:
#
#   julia --project=. -t auto examples/run_dayside_animation_1000.jl

using MarsASPEN

repo = normpath(joinpath(@__DIR__, ".."))
output = length(ARGS) >= 1 ? ARGS[1] : joinpath(
    repo, "examples", "output", "dayside_particle_animation_1000.mat",
)

# Use the same representative atmosphere selected by the other dayside
# examples. Changing these two keywords allows the animation to be repeated
# for another season or solar-activity condition.
model = load_model(repo; solar="solar_min", ls=0)

config = MonteCarloConfig(
    n_particles=1_000,
    injection_geometry=:dayside_uniform,
    initial_altitude_km=600.0,
    initial_speed_m_s=400_000.0,
    initial_charge_state=1,
    initial_temperature_ev=10.0,
    min_altitude_km=80.0,
    max_altitude_km=600.0,
    max_collisions=nothing,
    max_step_m=1_000.0,
    seed=73,
)

# Wn is a density weight in m^-3. It contains no velocity factor. The
# animation itself shows particle locations and charge states, while the
# weight is retained in the MAT metadata for reproducibility.
weighting = MonteCarloWeight(
    sampling_temperature_factor=1.0,
    source_number_density_m3=5.0e6,
)

summaries, histories = run_detailed_ensemble(
    model, config; weighting=weighting,
)
write_detailed_mat(
    output, summaries, histories; config=config, weighting=weighting,
)

println("particles=$(config.n_particles)")
println("mean_collisions=$(sum(s.n_collisions for s in summaries) / config.n_particles)")
println("output=$(abspath(output))")

# Single 400 km/s H+ trajectory.
#
# Run from the repository root:
#   julia --project=. examples/run_hplus_trajectory.jl
#
# This diagnostic follows one particle only. It does not assign a physical
# source density or a Monte Carlo macro-particle weight. Collisions themselves
# remain stochastic because MarsASPEN is a particle collision model.
#
# H+ can become H ENA through charge exchange, and stripping can return neutral
# H to H+. The saved history therefore exposes charge-state transitions,
# reaction altitudes, collision energy losses, and scattering.

using MarsASPEN

# Resolve paths from this script rather than using machine-specific paths.
repo = normpath(joinpath(@__DIR__, ".."))
output = length(ARGS) >= 1 ? ARGS[1] :
    joinpath(repo, "examples", "output", "single_hplus_400kms.mat")

# Assemble atmosphere, hot O, H and H+ collision data, and scattering tables.
model = load_model(repo; solar="solar_min", ls=0)

# State one selects H+. Zero temperature makes the source monoenergetic.
# The seed makes the stochastic trajectory reproducible.
config = MonteCarloConfig(
    n_particles=1,
    initial_altitude_km=600.0,
    initial_speed_m_s=400_000.0,
    initial_charge_state=1,
    initial_temperature_ev=0.0,
    seed=21,
)
# No source density normalization is appropriate for this one-particle plot.
summaries, histories = run_detailed_ensemble(model, config)
# Store every propagation and collision event as aligned MAT columns.
write_detailed_mat(output, summaries, histories; config=config)
println("output=$(abspath(output))")

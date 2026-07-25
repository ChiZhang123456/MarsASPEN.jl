# Single 400 km/s neutral H ENA trajectory.
#
# Run from the repository root:
#   julia --project=. examples/run_h_ena_trajectory.jl
#
# This diagnostic follows one particle only. It does not assign a physical
# source density or a Monte Carlo macro-particle weight. Collisions themselves
# remain stochastic because MarsASPEN is a particle collision model.
#
# The detailed MAT output contains one row for every transport step and every
# collision event. Important columns include time, position, velocity,
# altitude, kinetic energy, charge state, target species, reaction code, and
# scattering angle. This event-level format is suitable for one trajectory,
# but not for a large ensemble because its file size grows rapidly.

using MarsASPEN

# Resolve paths from this script so the example is portable.
repo = normpath(joinpath(@__DIR__, ".."))
output = length(ARGS) >= 1 ? ARGS[1] :
    joinpath(repo, "examples", "output", "single_h_ena_400kms.mat")

# Assemble atmosphere, hot O, cross sections, energy losses, and scattering.
model = load_model(repo; solar="solar_min", ls=0)

# State zero is neutral H. Zero temperature disables thermal perturbations, so
# the initial velocity is exactly 400 km/s toward negative MSO X. The seed
# controls all subsequent stochastic collision decisions.
config = MonteCarloConfig(
    n_particles=1,
    initial_altitude_km=600.0,
    initial_speed_m_s=400_000.0,
    initial_charge_state=0,
    initial_temperature_ev=0.0,
    seed=21,
)
# No physical density weight is needed for one illustrative trajectory.
summaries, histories = run_detailed_ensemble(model, config)
# Flatten the complete event history into Python-readable MAT columns.
write_detailed_mat(output, summaries, histories; config=config)
println("output=$(abspath(output))")

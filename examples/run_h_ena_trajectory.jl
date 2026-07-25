# Single 400 km/s neutral H ENA trajectory.
#
# Run from the repository root:
#   julia --project=. examples/run_h_ena_trajectory.jl
#
# This diagnostic follows one particle only. It does not assign a physical
# source density or a Monte Carlo macro-particle weight. Collisions themselves
# remain stochastic because MarsASPEN is a particle collision model.

using MarsASPEN

repo = normpath(joinpath(@__DIR__, ".."))
output = length(ARGS) >= 1 ? ARGS[1] :
    joinpath(repo, "examples", "output", "single_h_ena_400kms.mat")

model = load_model(repo; solar="solar_min", ls=0)
config = MonteCarloConfig(
    n_particles=1,
    initial_altitude_km=600.0,
    initial_speed_m_s=400_000.0,
    initial_charge_state=0,
    initial_temperature_ev=0.0,
    seed=21,
)
# The default unit weighting is sufficient for a single diagnostic trajectory.
summaries, histories = run_detailed_ensemble(model, config)
write_detailed_mat(output, summaries, histories; config=config)
println("output=$(abspath(output))")

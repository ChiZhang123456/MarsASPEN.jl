# Single 400 km/s H+ trajectory.
#
# Run from the repository root:
#   julia --project=. examples/run_hplus_trajectory.jl

using MarsASPEN

repo = normpath(joinpath(@__DIR__, ".."))
output = length(ARGS) >= 1 ? ARGS[1] :
    joinpath(repo, "examples", "output", "single_hplus_400kms.mat")

model = load_model(repo; solar="solar_min", ls=0)
config = MonteCarloConfig(
    n_particles=1,
    initial_altitude_km=600.0,
    initial_speed_m_s=400_000.0,
    initial_charge_state=1,
    initial_temperature_ev=0.0,
    seed=21,
)
summaries, histories = run_detailed_ensemble(model, config)
write_detailed_mat(output, summaries, histories; config=config)
println("output=$(abspath(output))")

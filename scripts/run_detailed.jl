# Run a small full-history ensemble and save every transport and collision row.
# Detailed histories grow rapidly, so this mode is intended for trajectory
# inspection rather than million-particle production simulations.
using MarsASPEN
using Printf

repo = normpath(joinpath(@__DIR__, ".."))
n = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 10
output = length(ARGS) >= 2 ? ARGS[2] :
    joinpath(repo, "output", "aspen_julia_$(n)p_detailed.mat")

atmosphere_dir = get(ENV, "MARSASPEN_ATMOSPHERE_DIR", joinpath(repo, "data", "atmosphere"))
model = load_model(repo; atmosphere_data_dir=atmosphere_dir)
config = MonteCarloConfig(n_particles=n, seed=21, include_hot_o=true)
weighting = MonteCarloWeight()
elapsed = @elapsed summaries, histories = run_detailed_ensemble(
    model, config; weighting=weighting,
)
write_detailed_mat(
    output, summaries, histories; config=config, weighting=weighting,
)

@printf("n_particles=%d\n", n)
@printf("elapsed_s=%.6f\n", elapsed)
@printf("n_events=%d\n", sum(length, histories))
@printf("n_collisions=%d\n", sum(r.n_collisions for r in summaries))
@printf("output=%s\n", abspath(output))

# Performance and reproducibility benchmark.
# Usage: julia --project=. -t auto scripts/benchmark.jl [n_particles] [single]
# A small serial run compiles the kernel before the timed measurement.
using MarsASPEN
using Statistics
using Printf

repo = normpath(joinpath(@__DIR__, ".."))
n = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 1000
threaded = length(ARGS) < 2 || ARGS[2] != "single"
atmosphere_dir = get(ENV, "MARSASPEN_ATMOSPHERE_DIR", joinpath(repo, "data", "atmosphere"))
model = load_model(repo; atmosphere_data_dir=atmosphere_dir)
cfg = MonteCarloConfig(n_particles=n)
weighting = MonteCarloWeight()

# Compile with a small run before measuring steady-state transport.
run_ensemble(
    model, MonteCarloConfig(n_particles=2);
    weighting=weighting, threaded=false,
)
elapsed = @elapsed rows = run_ensemble(
    model, cfg; weighting=weighting, threaded=threaded,
)

@printf("implementation=julia\n")
@printf("n_particles=%d\n", n)
@printf("threads=%d\n", threaded ? Threads.nthreads() : 1)
@printf("elapsed_s=%.6f\n", elapsed)
@printf("particles_per_s=%.3f\n", n / elapsed)
@printf("final_energy_mean_ev=%.6f\n", mean(r.final_energy_ev for r in rows))
@printf("n_collisions_mean=%.6f\n", mean(r.n_collisions for r in rows))
@printf("n_steps_mean=%.6f\n", mean(r.n_steps for r in rows))
for code in UInt8(1):UInt8(5)
    @printf("stop_%d=%d\n", code, count(r -> r.stop_code == code, rows))
end

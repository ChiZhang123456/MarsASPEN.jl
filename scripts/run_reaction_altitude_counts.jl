using MarsASPEN
using Printf
using Statistics

repo = normpath(joinpath(@__DIR__, ".."))
n = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 1_000_000
output = length(ARGS) >= 2 ? ARGS[2] :
    joinpath(repo, "output", "reaction_altitude_counts_$(n)p.csv")
bin_width_km = length(ARGS) >= 3 ? parse(Float64, ARGS[3]) : 10.0

model = load_model(repo; solar="solar_min", ls=0)
config = MonteCarloConfig(n_particles=n)
edges = collect(80.0:bin_width_km:1000.0)

run_binned_ensemble(
    model, MonteCarloConfig(n_particles=2);
    altitude_edges_km=edges,
)
elapsed = @elapsed result = run_binned_ensemble(
    model, config; altitude_edges_km=edges,
)

mkpath(dirname(abspath(output)))
open(output, "w") do io
    println(io, "altitude_low_km,altitude_high_km,altitude_center_km,state_change,ionization,lya,elastic,chemical_total,all_total")
    for ibin in axes(result.reaction_counts, 1)
        counts = result.reaction_counts[ibin, :]
        chemical_total = sum(counts[1:3])
        all_total = sum(counts)
        center = (edges[ibin] + edges[ibin+1]) / 2
        println(
            io,
            "$(edges[ibin]),$(edges[ibin+1]),$center,$(counts[1]),$(counts[2])," *
            "$(counts[3]),$(counts[4]),$chemical_total,$all_total",
        )
    end
end

@printf("n_particles=%d\n", n)
@printf("threads=%d\n", Threads.nthreads())
@printf("elapsed_s=%.6f\n", elapsed)
@printf("particles_per_s=%.3f\n", n / elapsed)
@printf("state_change_count=%d\n", sum(result.reaction_counts[:,1]))
@printf("ionization_count=%d\n", sum(result.reaction_counts[:,2]))
@printf("lya_count=%d\n", sum(result.reaction_counts[:,3]))
@printf("elastic_count=%d\n", sum(result.reaction_counts[:,4]))
@printf("final_energy_mean_ev=%.6f\n", mean(r.final_energy_ev for r in result.summaries))
@printf("output=%s\n", abspath(output))

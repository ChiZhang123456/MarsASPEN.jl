# Low-memory altitude-energy histogram for neutral H ENA and H+.
# Arguments: particle count, MAT output, altitude-bin width, number of
# logarithmic energy bins, and macro-particle weight.
# Each segment contributes particle_weight * ds to its phase-space bin.
using MarsASPEN
using MAT
using Printf
using Statistics

repo = normpath(joinpath(@__DIR__, ".."))
n = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 1_000_000
output = length(ARGS) >= 2 ? ARGS[2] :
    joinpath(repo, "output", "phase_space_$(n)p.mat")
altitude_bin_km = length(ARGS) >= 3 ? parse(Float64, ARGS[3]) : 1.0
energy_bins = length(ARGS) >= 4 ? parse(Int, ARGS[4]) : 100
particle_weight = length(ARGS) >= 5 ? parse(Float64, ARGS[5]) : 1.0

model = load_model(repo; solar="solar_min", ls=0)
config = MonteCarloConfig(n_particles=n, particle_weight=particle_weight)
altitude_edges = collect(80.0:altitude_bin_km:600.0)
energy_edges = 10.0 .^ range(
    log10(config.min_energy_ev), log10(1000.0), length=energy_bins+1,
)

run_phase_space_ensemble(
    model, MonteCarloConfig(n_particles=2);
    altitude_edges_km=altitude_edges,
    energy_edges_ev=energy_edges,
)
elapsed = @elapsed result = run_phase_space_ensemble(
    model, config;
    altitude_edges_km=altitude_edges,
    energy_edges_ev=energy_edges,
)

mkpath(dirname(abspath(output)))
matwrite(output, Dict(
    "format_version" => "marsaspen_phase_space_v1",
    "n_particles" => n,
    "seed" => config.seed,
    "ls_deg" => 0,
    "f107" => 70,
    "initial_altitude_km" => config.initial_altitude_km,
    "initial_speed_m_s" => config.initial_speed_m_s,
    "minimum_altitude_km" => config.min_altitude_km,
    "particle_weight" => config.particle_weight,
    "altitude_edges_km" => altitude_edges,
    "energy_edges_ev" => energy_edges,
    "charge_state_names" => ["H_ENA", "Hplus"],
    "path_length_m" => result.path_length_m,
))

@printf("n_particles=%d\n", n)
@printf("threads=%d\n", Threads.nthreads())
@printf("elapsed_s=%.6f\n", elapsed)
@printf("particles_per_s=%.3f\n", n / elapsed)
@printf("particle_weight=%.9g\n", config.particle_weight)
@printf("H_ENA_weighted_path_length_km=%.9g\n", sum(result.path_length_m[:,:,1]) / 1000)
@printf("Hplus_weighted_path_length_km=%.9g\n", sum(result.path_length_m[:,:,2]) / 1000)
@printf("final_energy_mean_ev=%.6f\n", mean(r.final_energy_ev for r in result.summaries))
@printf("output=%s\n", abspath(output))

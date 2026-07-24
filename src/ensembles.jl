function run_ensemble(model::AspenModel, cfg::MonteCarloConfig; threaded::Bool=true)
    out = Vector{ParticleSummary}(undef, cfg.n_particles)
    if threaded && nthreads() > 1
        @threads for i in eachindex(out)
            out[i] = first(run_particle_core(model, cfg, i, false))
        end
    else
        for i in eachindex(out)
            out[i] = first(run_particle_core(model, cfg, i, false))
        end
    end
    out
end
function run_binned_ensemble(
    model::AspenModel,
    cfg::MonteCarloConfig;
    altitude_edges_km::AbstractVector{<:Real}=collect(80.0:10.0:1000.0),
)
    edges = Float64.(altitude_edges_km)
    length(edges) >= 2 || throw(ArgumentError("altitude_edges_km needs at least two edges"))
    all(diff(edges) .> 0) || throw(ArgumentError("altitude_edges_km must increase"))
    summaries = Vector{ParticleSummary}(undef, cfg.n_particles)
    thread_counts = [
        zeros(Int64, length(edges)-1, 4) for _ in 1:Base.Threads.maxthreadid()
    ]
    @threads for i in eachindex(summaries)
        tid = threadid()
        summaries[i] = first(run_particle_core(
            model, cfg, i, false;
            altitude_edges_km=edges,
            reaction_counts=thread_counts[tid],
        ))
    end
    counts = reduce(+, thread_counts)
    (
        summaries=summaries,
        altitude_edges_km=edges,
        reaction_names=REACTION_NAMES,
        reaction_counts=counts,
    )
end

function run_phase_space_ensemble(
    model::AspenModel,
    cfg::MonteCarloConfig;
    altitude_edges_km::AbstractVector{<:Real}=collect(80.0:1.0:600.0),
    energy_edges_ev::AbstractVector{<:Real}=10.0 .^ range(1.0, 3.0, length=101),
)
    altitude_edges = Float64.(altitude_edges_km)
    energy_edges = Float64.(energy_edges_ev)
    length(altitude_edges) >= 2 ||
        throw(ArgumentError("altitude_edges_km needs at least two edges"))
    length(energy_edges) >= 2 ||
        throw(ArgumentError("energy_edges_ev needs at least two edges"))
    all(diff(altitude_edges) .> 0) ||
        throw(ArgumentError("altitude_edges_km must increase"))
    all(diff(energy_edges) .> 0) ||
        throw(ArgumentError("energy_edges_ev must increase"))

    summaries = Vector{ParticleSummary}(undef, cfg.n_particles)
    thread_weights = [
        zeros(
            Float64, length(altitude_edges)-1, length(energy_edges)-1, 2,
        ) for _ in 1:Base.Threads.maxthreadid()
    ]
    @threads for i in eachindex(summaries)
        tid = threadid()
        summaries[i] = first(run_particle_core(
            model, cfg, i, false;
            altitude_edges_km=altitude_edges,
            energy_edges_ev=energy_edges,
            path_length_m=thread_weights[tid],
        ))
    end
    (
        summaries=summaries,
        altitude_edges_km=altitude_edges,
        energy_edges_ev=energy_edges,
        charge_state_names=("H_ENA", "Hplus"),
        path_length_m=reduce(+, thread_weights),
    )
end

function run_detailed_ensemble(model::AspenModel, cfg::MonteCarloConfig)
    summaries = Vector{ParticleSummary}(undef, cfg.n_particles)
    histories = Vector{Vector{HistoryEvent}}(undef, cfg.n_particles)
    @threads for i in 1:cfg.n_particles
        summaries[i], histories[i] = run_particle_core(model, cfg, i, true)
    end
    summaries, histories
end

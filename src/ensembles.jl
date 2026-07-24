# Threaded ensemble drivers and low-memory diagnostics.
#
# Each particle has an independent RNG, so @threads scheduling does not change
# numerical results. Histogram modes allocate one accumulator per Julia thread
# and reduce them only after transport, avoiding locks in the inner loop.

"""Run compact particle summaries without storing trajectory histories."""
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
"""
Run an ensemble and count collision events in altitude bins.

The output matrix has dimensions `(altitude bin, reaction)`. Counts are raw
events, not physical rates, and one particle can contribute many events.
"""
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

"""
Accumulate H ENA and H+ altitude-energy distributions.

The returned array has dimensions `(altitude bin, energy bin, charge state)`.
Each segment contributes its particle density weight times `ds`. Without a
physical source density this is `particle_weight * ds`. This is a track-length
estimator, not an event-count histogram.
"""
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
    total_importance_weight = initial_importance_weight_sum(cfg)
    @threads for i in eachindex(summaries)
        tid = threadid()
        summaries[i] = first(run_particle_core(
            model, cfg, i, false;
            altitude_edges_km=altitude_edges,
            energy_edges_ev=energy_edges,
            path_length_m=thread_weights[tid],
            total_importance_weight=total_importance_weight,
        ))
    end
    (
        summaries=summaries,
        altitude_edges_km=altitude_edges,
        energy_edges_ev=energy_edges,
        charge_state_names=("H_ENA", "Hplus"),
        path_length_m=reduce(+, thread_weights),
        total_importance_weight=total_importance_weight,
    )
end

"""Compute the realized normalization `sum(f/fs)` for the source ensemble."""
function initial_importance_weight_sum(cfg::MonteCarloConfig)
    if cfg.initial_temperature_ev <= 0 || cfg.sampling_temperature_factor == 1
        return Float64(cfg.n_particles)
    end
    thread_sums = zeros(Float64, Base.Threads.maxthreadid())
    @threads for particle_id in 1:cfg.n_particles
        rng = Xoshiro(hash((cfg.seed, particle_id)))
        importance_weight = last(sample_initial_velocity(cfg, rng))
        thread_sums[threadid()] += importance_weight
    end
    sum(thread_sums)
end

"""
Accumulate upward and downward differential-flux numerators.

`flux[altitude surface, energy bin, charge, direction]` contains the sum of
injection flux weights for every crossing. Charge order is H ENA then H+.
Direction order is downward then upward. Divide by energy-bin width to obtain
`m^-2 s^-1 eV^-1` when physical density weighting is enabled.

This production driver intentionally does not retain per-particle summaries,
which keeps memory bounded for ensembles of ten million particles.
"""
function run_directional_flux_ensemble(
    model::AspenModel,
    cfg::MonteCarloConfig;
    altitude_surfaces_km::AbstractVector{<:Real}=collect(100.0:0.5:300.0),
    energy_edges_ev::AbstractVector{<:Real}=collect(10.0:2.0:1500.0),
)
    altitude_surfaces = Float64.(altitude_surfaces_km)
    energy_edges = Float64.(energy_edges_ev)
    all(diff(altitude_surfaces) .> 0) ||
        throw(ArgumentError("altitude_surfaces_km must increase"))
    all(diff(energy_edges) .> 0) ||
        throw(ArgumentError("energy_edges_ev must increase"))
    length(altitude_surfaces) >= 1 ||
        throw(ArgumentError("at least one altitude surface is required"))
    length(energy_edges) >= 2 ||
        throw(ArgumentError("energy_edges_ev needs at least two edges"))

    nslots = Base.Threads.maxthreadid()
    thread_flux = [
        zeros(
            Float64, length(altitude_surfaces), length(energy_edges)-1, 2, 2,
        ) for _ in 1:nslots
    ]
    thread_stop_counts = [zeros(Int64, 5) for _ in 1:nslots]
    thread_final_energy = zeros(Float64, nslots)
    thread_collisions = zeros(Int64, nslots)
    thread_steps = zeros(Int64, nslots)
    total_importance_weight = initial_importance_weight_sum(cfg)

    @threads for particle_id in 1:cfg.n_particles
        tid = threadid()
        summary = first(run_particle_core(
            model, cfg, particle_id, false;
            flux_altitude_km=altitude_surfaces,
            flux_energy_edges_ev=energy_edges,
            directional_flux=thread_flux[tid],
            total_importance_weight=total_importance_weight,
        ))
        1 <= summary.stop_code <= 5 &&
            (thread_stop_counts[tid][summary.stop_code] += 1)
        thread_final_energy[tid] += summary.final_energy_ev
        thread_collisions[tid] += summary.n_collisions
        thread_steps[tid] += summary.n_steps
    end

    (
        altitude_surfaces_km=altitude_surfaces,
        energy_edges_ev=energy_edges,
        charge_state_names=("H_ENA", "Hplus"),
        direction_names=("downward", "upward"),
        flux_m2_s=reduce(+, thread_flux),
        stop_counts=reduce(+, thread_stop_counts),
        final_energy_mean_ev=sum(thread_final_energy) / cfg.n_particles,
        collision_mean=sum(thread_collisions) / cfg.n_particles,
        step_mean=sum(thread_steps) / cfg.n_particles,
        total_importance_weight=total_importance_weight,
    )
end

"""Run full-history trajectories for small diagnostic ensembles."""
function run_detailed_ensemble(model::AspenModel, cfg::MonteCarloConfig)
    summaries = Vector{ParticleSummary}(undef, cfg.n_particles)
    histories = Vector{Vector{HistoryEvent}}(undef, cfg.n_particles)
    @threads for i in 1:cfg.n_particles
        summaries[i], histories[i] = run_particle_core(model, cfg, i, true)
    end
    summaries, histories
end

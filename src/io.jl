# Detailed trajectory serialization.
#
# Histories are flattened into column arrays because this layout is efficient
# in MAT/HDF5 and easy to read from Python. Particle start indices and event
# counts reconstruct each variable-length trajectory without object arrays.

"""
Write detailed particle histories to a compressed MAT v7.3 file.

Event rows use event type 0 initial, 1 transport step, 2 collision, and 3 final.
Collision rows include the state before and after the event, target code,
reaction code, energy loss, position, velocity, and charge state.
"""
function write_detailed_mat(filename::AbstractString, summaries, histories;
                            config::MonteCarloConfig,
                            weighting::MonteCarloWeight=MonteCarloWeight())
    # Flatten variable-length per-particle vectors into a table-like layout.
    events = reduce(vcat, histories; init=HistoryEvent[])
    starts = cumsum(vcat(1, length.(histories)[1:end-1]))
    data = Dict{String,Any}(
        "format_version" => "aspen_julia_history_v1",
        "event_type_names" => ["initial", "transport", "collision", "final"],
        "target_names" => ["none", "CO2", "O", "N2"],
        "reaction_names" => ["none", "state_change", "ionization", "lya", "elastic"],
        "particle_start_index_1based" => starts,
        "particle_event_count" => length.(histories),
        "particle_id" => getfield.(events, :particle_id),
        "event_type" => getfield.(events, :event_type),
        "step" => getfield.(events, :step),
        "collision_number" => getfield.(events, :collision),
        "time_s" => getfield.(events, :time_s),
        "x_m" => getfield.(events, :x_m), "y_m" => getfield.(events, :y_m),
        "z_m" => getfield.(events, :z_m),
        "vx_m_s" => getfield.(events, :vx_m_s), "vy_m_s" => getfield.(events, :vy_m_s),
        "vz_m_s" => getfield.(events, :vz_m_s),
        "altitude_km" => getfield.(events, :altitude_km),
        "energy_before_ev" => getfield.(events, :energy_before_ev),
        "energy_ev" => getfield.(events, :energy_ev),
        "vx_before_m_s" => getfield.(events, :vx_before_m_s),
        "vy_before_m_s" => getfield.(events, :vy_before_m_s),
        "vz_before_m_s" => getfield.(events, :vz_before_m_s),
        "charge_state" => getfield.(events, :charge_state),
        "target_code" => getfield.(events, :target),
        "reaction_code" => getfield.(events, :reaction),
        "energy_loss_ev" => getfield.(events, :energy_loss_ev),
        "final_energy_ev" => getfield.(summaries, :final_energy_ev),
        "final_altitude_km" => getfield.(summaries, :final_altitude_km),
        "n_steps" => getfield.(summaries, :n_steps),
        "n_collisions" => getfield.(summaries, :n_collisions),
        "stop_code" => getfield.(summaries, :stop_code),
        "config_n_particles" => config.n_particles,
        "config_seed" => config.seed,
        "config_include_hot_o" => config.include_hot_o,
        "config_initial_altitude_km" => config.initial_altitude_km,
        "config_initial_speed_m_s" => config.initial_speed_m_s,
        "config_initial_charge_state" => config.initial_charge_state,
        "config_initial_temperature_ev" => config.initial_temperature_ev,
        "weight_sampling_temperature_factor" =>
            weighting.sampling_temperature_factor,
        "weight_source_number_density_m3" =>
            weighting.source_number_density_m3,
        "weight_unit_particle_weight" => weighting.unit_particle_weight,
        "config_max_step_m" => config.max_step_m,
    )
    mkpath(dirname(abspath(filename)))
    matwrite(filename, data; compress=true)
    filename
end

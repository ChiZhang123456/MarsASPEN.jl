# Threaded ensemble drivers and low-memory diagnostics.
#
# Each particle has an independent RNG, so @threads scheduling does not change
# numerical results. Histogram modes allocate one accumulator per Julia thread
# and reduce them only after transport, avoiding locks in the inner loop.

"""Run compact particle summaries without storing trajectory histories."""
function run_ensemble(
    model::AspenModel,
    cfg::MonteCarloConfig;
    weighting::MonteCarloWeight=MonteCarloWeight(),
    threaded::Bool=true,
)
    out = Vector{ParticleSummary}(undef, cfg.n_particles)
    if threaded && nthreads() > 1
        @threads for i in eachindex(out)
            out[i] = first(run_particle_core(
                model, cfg, i, false; weighting=weighting,
            ))
        end
    else
        for i in eachindex(out)
            out[i] = first(run_particle_core(
                model, cfg, i, false; weighting=weighting,
            ))
        end
    end
    out
end

"""
Sample source positions, velocities, and weights without transporting them.

This diagnostic entry point uses exactly the same deterministic per-particle
random stream and source routines as `run_particle_core`. For a dayside
surface, the returned radial velocity is calculated separately at every
sampled position. The physical inward flux weight is `Wn * max(-Vr, 0)`.
"""
function sample_injection_ensemble(
    cfg::MonteCarloConfig;
    weighting::MonteCarloWeight=MonteCarloWeight(),
)
    n = cfg.n_particles
    position = zeros(Float64, n, 3)
    velocity = zeros(Float64, n, 3)
    longitude = zeros(Float64, n)
    latitude = zeros(Float64, n)
    sza = zeros(Float64, n)
    radial_velocity = zeros(Float64, n)
    density_weight = zeros(Float64, n)
    flux_weight = zeros(Float64, n)
    importance_weight = zeros(Float64, n)
    total_importance_weight = initial_importance_weight_sum(cfg, weighting)
    @threads for particle_id in 1:n
        rng = Xoshiro(hash((cfg.seed, particle_id)))
        charge, vx, vy, vz, wi =
            sample_initial_velocity(cfg, weighting, rng)
        x, y, z = sample_initial_position(cfg, rng)
        radius = sqrt(x*x + y*y + z*z)
        vr = (vx*x + vy*y + vz*z) / radius
        wn = weighting.source_number_density_m3 > 0 ?
            particle_density_weight(
                wi, weighting.source_number_density_m3,
                total_importance_weight,
            ) : weighting.unit_particle_weight
        coordinates = cartesian_to_lon_lat_alt(x, y, z)
        position[particle_id, :] .= (x, y, z)
        velocity[particle_id, :] .= (vx, vy, vz)
        longitude[particle_id] =
            mod(coordinates.lon_deg + 180.0, 360.0) - 180.0
        latitude[particle_id] = coordinates.lat_deg
        sza[particle_id] = rad2deg(acos(clamp(x / radius, -1.0, 1.0)))
        radial_velocity[particle_id] = vr
        density_weight[particle_id] = wn
        flux_weight[particle_id] = wn * max(-vr, 0.0)
        importance_weight[particle_id] = wi
        charge == cfg.initial_charge_state ||
            error("sampled charge state changed during injection")
    end
    (
        position_m=position,
        velocity_m_s=velocity,
        longitude_deg=longitude,
        latitude_deg=latitude,
        solar_zenith_angle_deg=sza,
        radial_velocity_m_s=radial_velocity,
        density_weight_m3=density_weight,
        inward_flux_weight_m2_s=flux_weight,
        importance_weight=importance_weight,
        total_importance_weight=total_importance_weight,
        injection_geometry=String(cfg.injection_geometry),
    )
end

"""
Run a dayside ensemble and accumulate complete three-dimensional diagnostics.

The macro-particle rate is obtained from the local inward injection flux and
the full dayside injection area. Within each spatial cell, residence time
gives number density, path length gives total scalar flux, and radial
displacement gives signed, upward, and downward radial flux. Realized
collisions give physical reaction rates resolved by charge, target, and
reaction channel.
"""
function run_spatial_grid_ensemble(
    model::AspenModel,
    cfg::MonteCarloConfig;
    weighting::MonteCarloWeight,
    altitude_edges_km::AbstractVector{<:Real}=collect(80.0:1.0:600.0),
)
    cfg.injection_geometry === :dayside_uniform ||
        throw(ArgumentError(
            "3D spatial grids require injection_geometry=:dayside_uniform",
        ))
    weighting.source_number_density_m3 > 0 ||
        throw(ArgumentError(
            "3D spatial grids require a positive source_number_density_m3",
        ))
    minimum(diff(Float64.(altitude_edges_km))) >= cfg.max_step_m / 1000 ||
        throw(ArgumentError(
            "max_step_m must not exceed the smallest altitude-bin width",
        ))
    grid = create_spatial_grid(
        model; altitude_edges_km=altitude_edges_km,
    )
    nslots = Base.Threads.maxthreadid()
    thread_stop_counts = [zeros(Int64, 5) for _ in 1:nslots]
    thread_collisions = zeros(Int64, nslots)
    thread_steps = zeros(Int64, nslots)
    total_importance_weight = initial_importance_weight_sum(cfg, weighting)
    @threads for particle_id in 1:cfg.n_particles
        tid = threadid()
        summary = first(run_particle_core(
            model, cfg, particle_id, false;
            spatial_grid=grid,
            weighting=weighting,
            total_importance_weight=total_importance_weight,
        ))
        1 <= summary.stop_code <= 5 &&
            (thread_stop_counts[tid][summary.stop_code] += 1)
        thread_collisions[tid] += summary.n_collisions
        thread_steps[tid] += summary.n_steps
    end
    injection_radius_m = (MARS_RADIUS_KM + cfg.initial_altitude_km) * 1000
    (
        grid=grid,
        stop_counts=reduce(+, thread_stop_counts),
        total_collisions=sum(thread_collisions),
        total_steps=sum(thread_steps),
        total_importance_weight=total_importance_weight,
        dayside_injection_area_m2=2pi * injection_radius_m^2,
    )
end
"""
Run an ensemble and count collision events in altitude bins.

The output matrix has dimensions `(altitude bin, reaction)`. Counts are raw
events, not physical rates, and one particle can contribute many events.
"""
function run_binned_ensemble(
    model::AspenModel,
    cfg::MonteCarloConfig;
    weighting::MonteCarloWeight=MonteCarloWeight(),
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
            weighting=weighting,
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
physical source density this is `unit_particle_weight * ds`. This is a track-length
estimator, not an event-count histogram.
"""
function run_phase_space_ensemble(
    model::AspenModel,
    cfg::MonteCarloConfig;
    weighting::MonteCarloWeight=MonteCarloWeight(),
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
    total_importance_weight = initial_importance_weight_sum(cfg, weighting)
    @threads for i in eachindex(summaries)
        tid = threadid()
        summaries[i] = first(run_particle_core(
            model, cfg, i, false;
            altitude_edges_km=altitude_edges,
            energy_edges_ev=energy_edges,
            path_length_m=thread_weights[tid],
            weighting=weighting,
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
function initial_importance_weight_sum(
    cfg::MonteCarloConfig, weighting::MonteCarloWeight,
)
    if cfg.initial_temperature_ev <= 0 ||
       weighting.sampling_temperature_factor == 1
        return Float64(cfg.n_particles)
    end
    thread_sums = zeros(Float64, Base.Threads.maxthreadid())
    @threads for particle_id in 1:cfg.n_particles
        rng = Xoshiro(hash((cfg.seed, particle_id)))
        importance_weight = last(sample_initial_velocity(cfg, weighting, rng))
        thread_sums[threadid()] += importance_weight
    end
    sum(thread_sums)
end

"""
Sum density weights when particles cross fixed altitude surfaces.

The returned array has dimensions `(altitude surface, energy bin, charge)`.
Every upward or downward crossing contributes the particle density weight
exactly once. No path-length factor and no division by energy-bin width are
applied. Charge order is H ENA then H+.
"""
function run_density_crossing_ensemble(
    model::AspenModel,
    cfg::MonteCarloConfig;
    weighting::MonteCarloWeight=MonteCarloWeight(),
    altitude_surfaces_km::AbstractVector{<:Real}=collect(80.5:1.0:599.5),
    energy_edges_ev::AbstractVector{<:Real}=10.0 .^ range(0.0, 4.0, length=101),
)
    altitude_surfaces = Float64.(altitude_surfaces_km)
    energy_edges = Float64.(energy_edges_ev)
    all(diff(altitude_surfaces) .> 0) ||
        throw(ArgumentError("altitude_surfaces_km must increase"))
    all(diff(energy_edges) .> 0) ||
        throw(ArgumentError("energy_edges_ev must increase"))

    thread_density = [
        zeros(Float64, length(altitude_surfaces), length(energy_edges)-1, 2)
        for _ in 1:Base.Threads.maxthreadid()
    ]
    total_importance_weight = initial_importance_weight_sum(cfg, weighting)
    @threads for particle_id in 1:cfg.n_particles
        tid = threadid()
        run_particle_core(
            model, cfg, particle_id, false;
            flux_altitude_km=altitude_surfaces,
            flux_energy_edges_ev=energy_edges,
            density_crossings=thread_density[tid],
            weighting=weighting,
            total_importance_weight=total_importance_weight,
        )
    end
    (
        altitude_surfaces_km=altitude_surfaces,
        energy_edges_ev=energy_edges,
        charge_state_names=("H_ENA", "Hplus"),
        density_weight_sum_m3=reduce(+, thread_density),
        total_importance_weight=total_importance_weight,
    )
end

"""
Accumulate local radial number flux on spherical altitude surfaces.

For every surface crossing, the contribution is `Wn * abs(Vr)`, where `Wn`
has units m^-3 and `Vr = v dot r_hat` has units m s^-1. The result therefore
has units m^-2 s^-1. Array dimensions are `(altitude, charge, direction)`.
Charge order is H ENA then H+, and direction order is downward then upward.

The signed outward radial flux is `upward - downward`. The net downward flux
is the opposite sign, `downward - upward`.
"""
function run_radial_flux_ensemble(
    model::AspenModel,
    cfg::MonteCarloConfig;
    weighting::MonteCarloWeight=MonteCarloWeight(),
    altitude_surfaces_km::AbstractVector{<:Real}=collect(80.5:1.0:599.5),
)
    altitude_surfaces = Float64.(altitude_surfaces_km)
    length(altitude_surfaces) >= 1 ||
        throw(ArgumentError("at least one altitude surface is required"))
    all(diff(altitude_surfaces) .> 0) ||
        throw(ArgumentError("altitude_surfaces_km must increase"))

    thread_flux = [
        zeros(Float64, length(altitude_surfaces), 2, 2)
        for _ in 1:Base.Threads.maxthreadid()
    ]
    total_importance_weight = initial_importance_weight_sum(cfg, weighting)
    @threads for particle_id in 1:cfg.n_particles
        tid = threadid()
        run_particle_core(
            model, cfg, particle_id, false;
            flux_altitude_km=altitude_surfaces,
            radial_flux=thread_flux[tid],
            weighting=weighting,
            total_importance_weight=total_importance_weight,
        )
    end
    flux = reduce(+, thread_flux)
    (
        altitude_surfaces_km=altitude_surfaces,
        charge_state_names=("H_ENA", "Hplus"),
        direction_names=("downward", "upward"),
        radial_flux_m2_s=flux,
        signed_outward_flux_m2_s=flux[:, :, 2] .- flux[:, :, 1],
        net_downward_flux_m2_s=flux[:, :, 1] .- flux[:, :, 2],
        total_importance_weight=total_importance_weight,
    )
end

"""
Estimate target ionization rates from local radial crossing fluxes.

At each spherical altitude crossing, the transport kernel accumulates

`Wn * abs(Vr) * sigma_ion(E)`

separately by projectile charge state and direction. `Wn` is in m^-3, `Vr`
is in m s^-1, and the ionization cross section is in m^2, so this intermediate
coefficient is in s^-1. The driver multiplies it by the selected neutral
target density in m^-3 to obtain a volume ionization rate in m^-3 s^-1.

The target density is evaluated at `density_lon_deg` and `density_lat_deg`.
For target O, `cfg.include_hot_o` controls whether MAMPS hot O is added to
MGITM cold O. Charge order is H ENA then H+, and direction order is downward
then upward.
"""
function run_target_ionization_rate_ensemble(
    model::AspenModel,
    cfg::MonteCarloConfig;
    weighting::MonteCarloWeight=MonteCarloWeight(),
    altitude_surfaces_km::AbstractVector{<:Real}=collect(80.5:1.0:599.5),
    target::Symbol=:O,
    density_lon_deg::Real=0.0,
    density_lat_deg::Real=0.0,
)
    target_index = target === :CO2 ? 1 :
                   target === :O ? 2 :
                   target === :N2 ? 3 :
                   throw(ArgumentError("target must be :CO2, :O, or :N2"))
    altitude_surfaces = Float64.(altitude_surfaces_km)
    length(altitude_surfaces) >= 1 ||
        throw(ArgumentError("at least one altitude surface is required"))
    all(diff(altitude_surfaces) .> 0) ||
        throw(ArgumentError("altitude_surfaces_km must increase"))

    thread_flux_sigma = [
        zeros(Float64, length(altitude_surfaces), 2, 2)
        for _ in 1:Base.Threads.maxthreadid()
    ]
    thread_radial_flux = [
        zeros(Float64, length(altitude_surfaces), 2, 2)
        for _ in 1:Base.Threads.maxthreadid()
    ]
    total_importance_weight = initial_importance_weight_sum(cfg, weighting)
    @threads for particle_id in 1:cfg.n_particles
        tid = threadid()
        run_particle_core(
            model, cfg, particle_id, false;
            flux_altitude_km=altitude_surfaces,
            radial_flux=thread_radial_flux[tid],
            ionization_flux_sigma=thread_flux_sigma[tid],
            ionization_target=target_index,
            weighting=weighting,
            total_importance_weight=total_importance_weight,
        )
    end
    flux_sigma = reduce(+, thread_flux_sigma)
    radial_flux = reduce(+, thread_radial_flux)

    target_density = Vector{Float64}(undef, length(altitude_surfaces))
    for (ia, altitude_km) in pairs(altitude_surfaces)
        density = neutral_density(
            model, density_lon_deg, density_lat_deg, altitude_km;
            include_hot_o=cfg.include_hot_o,
        )
        target_density[ia] = target === :CO2 ? density.CO2 :
                             target === :O ? density.O : density.N2
    end
    rate = flux_sigma .* reshape(target_density, :, 1, 1)
    (
        altitude_surfaces_km=altitude_surfaces,
        target_name=String(target),
        target_density_m3=target_density,
        density_lon_deg=Float64(density_lon_deg),
        density_lat_deg=Float64(density_lat_deg),
        charge_state_names=("H_ENA", "Hplus"),
        direction_names=("downward", "upward"),
        radial_flux_m2_s=radial_flux,
        flux_times_ionization_cross_section_s1=flux_sigma,
        ionization_rate_m3_s1=rate,
        ionization_rate_by_charge_m3_s1=dropdims(sum(rate; dims=3); dims=3),
        total_ionization_rate_m3_s1=vec(sum(rate; dims=(2, 3))),
        total_importance_weight=total_importance_weight,
    )
end

"""Mean fractional projectile-energy loss in one sampled elastic collision."""
function mean_elastic_loss_fractions(model::AspenModel)
    fractions = zeros(Float64, 2, 3)
    random_grid = model.cross_sections.scatter_r
    angles = deg2rad.(model.cross_sections.scatter_theta)
    for charge in 0:1, target in 1:3
        projectile_mass = charge == 1 ? HP_MASS : H_MASS
        ratio = projectile_mass / TARGET_MASS[target]
        loss = similar(angles)
        for i in eachindex(angles)
            theta = angles[i]
            disc = max(1 - (ratio * sin(theta))^2, 0.0)
            speed_ratio = max(
                (ratio * cos(theta) + sqrt(disc)) / (1 + ratio), 0.0,
            )
            loss[i] = 1 - speed_ratio^2
        end
        integral = sum(
            0.5 * (loss[i] + loss[i + 1]) *
            (random_grid[i + 1] - random_grid[i])
            for i in 1:length(random_grid)-1
        )
        fractions[charge + 1, target] =
            integral / (random_grid[end] - random_grid[1])
    end
    fractions
end

"""
Calculate point-source diagnostic maps on one spherical altitude surface.

The first four output fields use the local total projectile speed, not radial
speed. Array order is longitude bin, latitude bin, diagnostic. Diagnostics are
O ionization, CO2 ionization, Ly-alpha volume emission, and projectile
collision-energy transfer. The terminal-energy field adds the remaining
energy of particles that stop below `cfg.min_energy_ev` inside the requested
altitude shell.

Because the present source starts at one geographic point, these arrays are a
local footprint rather than a globally normalized dayside map.
"""
function run_surface_diagnostic_map_ensemble(
    model::AspenModel,
    cfg::MonteCarloConfig;
    weighting::MonteCarloWeight=MonteCarloWeight(),
    altitude_km::Real=120.0,
    longitude_edges_deg::AbstractVector{<:Real}=collect(-30.0:1.0:30.0),
    latitude_edges_deg::AbstractVector{<:Real}=collect(-30.0:1.0:30.0),
    thermalization_shell_half_width_km::Real=0.5,
)
    lon_edges = Float64.(longitude_edges_deg)
    lat_edges = Float64.(latitude_edges_deg)
    all(diff(lon_edges) .> 0) ||
        throw(ArgumentError("longitude_edges_deg must increase"))
    all(diff(lat_edges) .> 0) ||
        throw(ArgumentError("latitude_edges_deg must increase"))
    half_width = Float64(thermalization_shell_half_width_km)
    half_width > 0 ||
        throw(ArgumentError("thermalization_shell_half_width_km must be positive"))
    surface_altitude = Float64(altitude_km)
    shell_edges = (
        surface_altitude - half_width, surface_altitude + half_width,
    )
    thread_maps = [
        zeros(Float64, length(lon_edges)-1, length(lat_edges)-1, 4)
        for _ in 1:Base.Threads.maxthreadid()
    ]
    thread_thermal = [
        zeros(Float64, length(lon_edges)-1, length(lat_edges)-1)
        for _ in 1:Base.Threads.maxthreadid()
    ]
    elastic_fractions = mean_elastic_loss_fractions(model)
    total_importance_weight = initial_importance_weight_sum(cfg, weighting)
    @threads for particle_id in 1:cfg.n_particles
        tid = threadid()
        run_particle_core(
            model, cfg, particle_id, false;
            flux_altitude_km=[surface_altitude],
            surface_diagnostics=thread_maps[tid],
            surface_altitude_km=surface_altitude,
            surface_longitude_edges_deg=lon_edges,
            surface_latitude_edges_deg=lat_edges,
            elastic_loss_fraction=elastic_fractions,
            thermalization_energy_map=thread_thermal[tid],
            thermalization_shell_edges_km=shell_edges,
            weighting=weighting,
            total_importance_weight=total_importance_weight,
        )
    end
    diagnostics = reduce(+, thread_maps)
    thermalization = reduce(+, thread_thermal)
    collision_energy = diagnostics[:, :, 4]
    total_energy = collision_energy .+ thermalization
    (
        altitude_km=surface_altitude,
        longitude_edges_deg=lon_edges,
        latitude_edges_deg=lat_edges,
        diagnostic_names=(
            "O_ionization", "CO2_ionization", "Ly_alpha",
            "collision_energy_transfer",
        ),
        oxygen_ionization_rate_m3_s1=diagnostics[:, :, 1],
        co2_ionization_rate_m3_s1=diagnostics[:, :, 2],
        total_ionization_rate_m3_s1=
            diagnostics[:, :, 1] .+ diagnostics[:, :, 2],
        lya_volume_emission_rate_photons_m3_s1=diagnostics[:, :, 3],
        collision_energy_transfer_ev_m3_s1=collision_energy,
        thermalized_below_cutoff_ev_m3_s1=thermalization,
        total_energy_transfer_ev_m3_s1=total_energy,
        total_energy_transfer_w_m3=total_energy .* QE,
        elastic_mean_loss_fraction=elastic_fractions,
        thermalization_shell_edges_km=collect(shell_edges),
        total_importance_weight=total_importance_weight,
        source_geometry="single point at lon=0 deg, lat=0 deg, 600 km",
    )
end

"""
Estimate H Ly-alpha volume emission and radiative-energy profiles.

At every spherical altitude crossing, the transport kernel evaluates the
local total projectile speed and energy, then accumulates

`Wn * speed * sigma_Lya(E)`

by projectile charge state, target species, and direction. The effective
Ly-alpha production cross section is reaction channel 3. Multiplication by
the local CO2, O, or N2 density gives photons m^-3 s^-1, assuming the effective
cross section already includes the photon yield.

O density includes MAMPS hot O when `cfg.include_hot_o=true`. The radiative
energy rate multiplies the photon rate by `h*c/lambda` at 121.567 nm. It is a
radiative source term, not a local thermal-heating rate.
"""
function run_lya_volume_emission_ensemble(
    model::AspenModel,
    cfg::MonteCarloConfig;
    weighting::MonteCarloWeight=MonteCarloWeight(),
    altitude_surfaces_km::AbstractVector{<:Real}=collect(80.5:1.0:599.5),
    density_lon_deg::Real=0.0,
    density_lat_deg::Real=0.0,
)
    altitude_surfaces = Float64.(altitude_surfaces_km)
    length(altitude_surfaces) >= 1 ||
        throw(ArgumentError("at least one altitude surface is required"))
    all(diff(altitude_surfaces) .> 0) ||
        throw(ArgumentError("altitude_surfaces_km must increase"))

    thread_coefficients = [
        zeros(Float64, length(altitude_surfaces), 2, 3, 2)
        for _ in 1:Base.Threads.maxthreadid()
    ]
    total_importance_weight = initial_importance_weight_sum(cfg, weighting)
    @threads for particle_id in 1:cfg.n_particles
        tid = threadid()
        run_particle_core(
            model, cfg, particle_id, false;
            flux_altitude_km=altitude_surfaces,
            lya_rate_coefficient=thread_coefficients[tid],
            weighting=weighting,
            total_importance_weight=total_importance_weight,
        )
    end
    coefficient = reduce(+, thread_coefficients)

    target_density = Matrix{Float64}(undef, length(altitude_surfaces), 3)
    for (ia, altitude_km) in pairs(altitude_surfaces)
        density = neutral_density(
            model, density_lon_deg, density_lat_deg, altitude_km;
            include_hot_o=cfg.include_hot_o,
        )
        target_density[ia, :] .= (density.CO2, density.O, density.N2)
    end
    volume_emission = coefficient .* reshape(
        target_density, length(altitude_surfaces), 1, 3, 1,
    )
    by_charge = dropdims(sum(volume_emission; dims=(3, 4)); dims=(3, 4))
    by_target = dropdims(sum(volume_emission; dims=(2, 4)); dims=(2, 4))
    total_emission = vec(sum(volume_emission; dims=(2, 3, 4)))
    photon_energy_j = PLANCK_J_S * LIGHT_SPEED_M_S / LYA_WAVELENGTH_M
    (
        altitude_surfaces_km=altitude_surfaces,
        charge_state_names=("H_ENA", "Hplus"),
        target_names=("CO2", "O", "N2"),
        direction_names=("downward", "upward"),
        target_density_m3=target_density,
        density_lon_deg=Float64(density_lon_deg),
        density_lat_deg=Float64(density_lat_deg),
        lya_rate_coefficient_s1=coefficient,
        volume_emission_rate_photons_m3_s1=volume_emission,
        volume_emission_rate_by_charge_photons_m3_s1=by_charge,
        volume_emission_rate_by_target_photons_m3_s1=by_target,
        total_volume_emission_rate_photons_m3_s1=total_emission,
        photon_energy_j=photon_energy_j,
        radiative_energy_rate_w_m3=volume_emission .* photon_energy_j,
        radiative_energy_rate_by_charge_w_m3=by_charge .* photon_energy_j,
        total_radiative_energy_rate_w_m3=total_emission .* photon_energy_j,
        total_importance_weight=total_importance_weight,
    )
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
    weighting::MonteCarloWeight=MonteCarloWeight(),
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
    total_importance_weight = initial_importance_weight_sum(cfg, weighting)

    @threads for particle_id in 1:cfg.n_particles
        tid = threadid()
        summary = first(run_particle_core(
            model, cfg, particle_id, false;
            flux_altitude_km=altitude_surfaces,
            flux_energy_edges_ev=energy_edges,
            directional_flux=thread_flux[tid],
            weighting=weighting,
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
function run_detailed_ensemble(
    model::AspenModel,
    cfg::MonteCarloConfig;
    weighting::MonteCarloWeight=MonteCarloWeight(),
)
    summaries = Vector{ParticleSummary}(undef, cfg.n_particles)
    histories = Vector{Vector{HistoryEvent}}(undef, cfg.n_particles)
    @threads for i in 1:cfg.n_particles
        summaries[i], histories[i] = run_particle_core(
            model, cfg, i, true; weighting=weighting,
        )
    end
    summaries, histories
end

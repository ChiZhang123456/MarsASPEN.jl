# Single-particle Monte Carlo transport kernel.
#
# The algorithm samples a collision optical depth `-log(U)`, advances the
# projectile with adaptive spatial steps, and accumulates `alpha * ds` until
# the sampled optical depth is reached. A target and reaction are then drawn
# from local `n * sigma` weights. Energy, direction, and charge state are
# updated immediately before the next collision probability is evaluated.

"""
Rotate a velocity vector by polar angle `theta` and azimuth `phi`.

The original velocity defines the local polar axis. The function returns a
vector with magnitude `speed_after`, so scattering direction and energy loss
are handled independently.
"""
@inline function rotate_velocity(vx, vy, vz, speed_after, theta, phi)
    speed = sqrt(vx*vx + vy*vy + vz*vz)
    e0x, e0y, e0z = vx/speed, vy/speed, vz/speed
    hx, hy, hz = abs(e0z) > 0.9 ? (0.0, 1.0, 0.0) : (0.0, 0.0, 1.0)
    e1x, e1y, e1z = e0y*hz-e0z*hy, e0z*hx-e0x*hz, e0x*hy-e0y*hx
    en = sqrt(e1x*e1x + e1y*e1y + e1z*e1z)
    e1x, e1y, e1z = e1x/en, e1y/en, e1z/en
    e2x, e2y, e2z = e0y*e1z-e0z*e1y, e0z*e1x-e0x*e1z, e0x*e1y-e0y*e1x
    ct, st, cp, sp = cos(theta), sin(theta), cos(phi), sin(phi)
    speed_after * (ct*e0x + st*(cp*e1x + sp*e2x)),
    speed_after * (ct*e0y + st*(cp*e1y + sp*e2y)),
    speed_after * (ct*e0z + st*(cp*e1z + sp*e2z))
end

"""
Return position and signed radial velocity where a segment crosses a sphere.

The segment is `p(t) = p0 + t * (p1 - p0)`, with `0 <= t <= 1`. The velocity
is constant during one free-flight segment. Positive radial velocity is
outward or upward, and negative radial velocity is inward or downward.
"""
@inline function state_at_surface_crossing(
    x, y, z, xnew, ynew, znew, vx, vy, vz, altitude_km,
)
    dx, dy, dz = xnew - x, ynew - y, znew - z
    radius_m = (MARS_RADIUS_KM + altitude_km) * 1000
    a = dx*dx + dy*dy + dz*dz
    b = 2 * (x*dx + y*dy + z*dz)
    c = x*x + y*y + z*z - radius_m*radius_m
    discriminant = max(b*b - 4*a*c, 0.0)
    root = sqrt(discriminant)
    t1 = (-b - root) / (2a)
    t2 = (-b + root) / (2a)
    t = 0.0 <= t1 <= 1.0 ? t1 : clamp(t2, 0.0, 1.0)
    xc, yc, zc = x + t*dx, y + t*dy, z + t*dz
    vr = (vx*xc + vy*yc + vz*zc) / radius_m
    xc, yc, zc, vr
end

"""Return signed radial velocity at a spherical-surface crossing."""
@inline function radial_velocity_at_surface(
    x, y, z, xnew, ynew, znew, vx, vy, vz, altitude_km,
)
    _, _, _, vr = state_at_surface_crossing(
        x, y, z, xnew, ynew, znew, vx, vy, vz, altitude_km,
    )
    vr
end

"""
Sample the initial drifting-Maxwellian projectile velocity and importance weight.

The bulk velocity is along MSO -X. Thermal components are sampled at
`T_sample = sampling_temperature_factor * T`. The returned importance weight
is `f(T) / f(T_sample)`.
"""
@inline function sample_initial_velocity(
    cfg::MonteCarloConfig, weighting::MonteCarloWeight, rng,
)
    charge = cfg.initial_charge_state
    charge in (0, 1) || throw(ArgumentError("initial_charge_state must be 0 or 1"))
    mass = charge == 1 ? HP_MASS : H_MASS
    weighting.sampling_temperature_factor > 0 ||
        throw(ArgumentError("sampling_temperature_factor must be positive"))
    sampled_temperature_ev =
        max(cfg.initial_temperature_ev, 0.0) *
        weighting.sampling_temperature_factor
    thermal_sigma = sqrt(sampled_temperature_ev * QE / mass)
    if thermal_sigma > 0
        vx = -cfg.initial_speed_m_s + thermal_sigma * randn(rng)
        vy = thermal_sigma * randn(rng)
        vz = thermal_sigma * randn(rng)
    else
        # Preserve the original monoenergetic random stream when T=0.
        vx, vy, vz = -cfg.initial_speed_m_s, 0.0, 0.0
    end
    importance_weight = cfg.initial_temperature_ev > 0 ?
        maxwellian_importance_weight_3d(
            (-cfg.initial_speed_m_s, 0.0, 0.0), (vx, vy, vz),
            cfg.initial_temperature_ev, sampled_temperature_ev;
            mass_kg=mass,
        ) : 1.0
    charge, vx, vy, vz, importance_weight
end

"""
Return initial state, density weight, and physical crossing-flux weight.

When a source density is supplied, the density weight follows py_aspen:
`Wn_i = n_source * (f/fs)_i / sum(f/fs)`. Multiplying by the inward normal
speed gives the particle crossing-flux weight in m^-2 s^-1.
"""
@inline function initial_particle_state(
    cfg::MonteCarloConfig,
    weighting::MonteCarloWeight,
    rng,
    total_importance_weight::Float64,
)
    charge, vx, vy, vz, importance_weight =
        sample_initial_velocity(cfg, weighting, rng)
    density_weight = if weighting.source_number_density_m3 > 0
        density_weight = particle_density_weight(
            importance_weight, weighting.source_number_density_m3,
            total_importance_weight,
        )
        density_weight
    else
        weighting.unit_particle_weight
    end
    flux_weight = weighting.source_number_density_m3 > 0 ?
                  density_weight * max(-vx, 0.0) : density_weight
    charge, vx, vy, vz, density_weight, flux_weight
end
"""
Transport one particle until a physical or numerical stopping condition.

Optional accumulators support three output modes without duplicating the
physics kernel:

* `record=true` stores full step and collision history.
* `reaction_counts` accumulates collision events by altitude and reaction.
* `path_length_m` accumulates the particle density weight times `ds` by altitude, energy,
  and charge state.
* `directional_flux` counts weighted upward and downward crossings of altitude
  surfaces by energy and charge state.
* `density_crossings` sums particle density weights at altitude-surface
  crossings, without multiplying by path length or energy-bin width.
* `radial_flux` sums `Wn * abs(Vr)` at each crossing, separately for downward
  and upward directions. Its dimensions are altitude, charge, direction.
* `ionization_flux_sigma` sums `Wn * abs(Vr) * sigma_ion(E)` at crossings.
  Multiplication by the selected target density is performed by the ensemble
  driver to obtain an ionization rate in m^-3 s^-1.
* `lya_rate_coefficient` sums `Wn * speed * sigma_Lya(E)` by altitude,
  charge, target, and direction. Multiplication by each target density gives
  the Ly-alpha volume emission rate in photons m^-3 s^-1.

Particles begin at 600 km with the configured charge state and drifting
Maxwellian velocity. The random stream is derived from `(cfg.seed, id)`, which
makes results reproducible and independent of thread scheduling.
"""
function run_particle_core(
    model::AspenModel, cfg::MonteCarloConfig, id::Int, record::Bool;
    altitude_edges_km::Union{Nothing,Vector{Float64}}=nothing,
    reaction_counts::Union{Nothing,Matrix{Int64}}=nothing,
    energy_edges_ev::Union{Nothing,Vector{Float64}}=nothing,
    path_length_m::Union{Nothing,Array{Float64,3}}=nothing,
    flux_altitude_km::Union{Nothing,Vector{Float64}}=nothing,
    flux_energy_edges_ev::Union{Nothing,Vector{Float64}}=nothing,
    directional_flux::Union{Nothing,Array{Float64,4}}=nothing,
    density_crossings::Union{Nothing,Array{Float64,3}}=nothing,
    radial_flux::Union{Nothing,Array{Float64,3}}=nothing,
    ionization_flux_sigma::Union{Nothing,Array{Float64,3}}=nothing,
    ionization_target::Int=2,
    lya_rate_coefficient::Union{Nothing,Array{Float64,4}}=nothing,
    surface_diagnostics::Union{Nothing,Array{Float64,3}}=nothing,
    surface_altitude_km::Float64=120.0,
    surface_longitude_edges_deg::Union{Nothing,Vector{Float64}}=nothing,
    surface_latitude_edges_deg::Union{Nothing,Vector{Float64}}=nothing,
    elastic_loss_fraction::Union{Nothing,Matrix{Float64}}=nothing,
    thermalization_energy_map::Union{Nothing,Matrix{Float64}}=nothing,
    thermalization_shell_edges_km::NTuple{2,Float64}=(119.5, 120.5),
    weighting::MonteCarloWeight=MonteCarloWeight(),
    total_importance_weight::Float64=Float64(cfg.n_particles),
)
    # A separate deterministic RNG per particle prevents thread-order effects.
    rng = Xoshiro(hash((cfg.seed, id)))
    x, y, z = (MARS_RADIUS_KM + cfg.initial_altitude_km) * 1000, 0.0, 0.0
    charge, vx, vy, vz, particle_density_weight_m3, injection_flux_weight =
        initial_particle_state(cfg, weighting, rng, total_importance_weight)
    # `tau` is accumulated optical depth and
    # `threshold` is the exponentially distributed optical depth to collision.
    tau, threshold = 0.0, -log(rand(rng))
    steps = collisions = 0
    counts = zeros(Int, 4)
    stop = UInt8(0)
    elapsed_time = 0.0
    history = HistoryEvent[]
    if record
        alt0 = sqrt(x*x+y*y+z*z)/1000 - MARS_RADIUS_KM
        push!(history, HistoryEvent(id, 0, 0, 0, 0.0, x,y,z,vx,vy,vz,alt0,
            energy(vx,vy,vz,charge),energy(vx,vy,vz,charge),vx,vy,vz,
            Int8(charge), 0, 0, 0.0))
    end
    while collisions < cfg.max_collisions
        e = energy(vx, vy, vz, charge)
        alt, n, alpha = local_state(model, x, y, z, charge, e, cfg.include_hot_o)
        if e < cfg.min_energy_ev
            stop = 1; break
        elseif alt < cfg.min_altitude_km
            stop = 2; break
        elseif alt > cfg.max_altitude_km
            stop = 3; break
        elseif steps >= cfg.max_steps_per_collision * max(collisions + 1, 1)
            stop = 4; break
        end
        speed = sqrt(vx*vx + vy*vy + vz*vz)
        # Limit both optical-depth change and absolute spatial step length.
        candidate = alpha > 0 ? min(cfg.safety_factor / alpha, cfg.max_step_m) : cfg.max_step_m
        remaining = threshold - tau
        ds = alpha > 0 ? min(candidate, remaining / alpha) : candidate
        if !isnothing(path_length_m)
            # Midpoint binning reduces altitude error for a finite segment.
            # Weighting by ds avoids bias from adaptive step subdivision.
            mx = x + 0.5 * vx / speed * ds
            my = y + 0.5 * vy / speed * ds
            mz = z + 0.5 * vz / speed * ds
            midpoint_altitude = sqrt(mx*mx + my*my + mz*mz)/1000 - MARS_RADIUS_KM
            ia = searchsortedlast(altitude_edges_km, midpoint_altitude)
            ie = searchsortedlast(energy_edges_ev, e)
            if 1 <= ia < length(altitude_edges_km) &&
               1 <= ie < length(energy_edges_ev)
                path_length_m[ia, ie, charge + 1] +=
                    particle_density_weight_m3 * ds
            end
        end
        xnew = x + vx / speed * ds
        ynew = y + vy / speed * ds
        znew = z + vz / speed * ds
        if !isnothing(directional_flux) || !isnothing(density_crossings) ||
           !isnothing(radial_flux) || !isnothing(ionization_flux_sigma) ||
           !isnothing(lya_rate_coefficient) || !isnothing(surface_diagnostics)
            # Count crossings of fixed spherical altitude surfaces. Direction
            # 1 is downward and direction 2 is upward.
            altitude_after = sqrt(xnew*xnew + ynew*ynew + znew*znew)/1000 -
                             MARS_RADIUS_KM
            has_energy_accumulator =
                !isnothing(directional_flux) || !isnothing(density_crossings)
            ie = has_energy_accumulator ?
                searchsortedlast(flux_energy_edges_ev, e) : 0
            valid_energy = has_energy_accumulator &&
                1 <= ie < length(flux_energy_edges_ev)
            sigma_ion_segment = !isnothing(ionization_flux_sigma) ?
                sigma_at(
                    model.cross_sections, charge, ionization_target, 2, e,
                ) : 0.0
            sigma_lya_segment = !isnothing(lya_rate_coefficient) ?
                ntuple(
                    target -> sigma_at(
                        model.cross_sections, charge, target, 3, e,
                    ),
                    3,
                ) : (0.0, 0.0, 0.0)
            if altitude_after < alt
                first_surface = searchsortedlast(flux_altitude_km, altitude_after) + 1
                last_surface = searchsortedlast(flux_altitude_km, alt)
                for ia in first_surface:last_surface
                    if valid_energy
                        if !isnothing(directional_flux)
                            directional_flux[ia, ie, charge + 1, 1] +=
                                injection_flux_weight
                        end
                        if !isnothing(density_crossings)
                            density_crossings[ia, ie, charge + 1] +=
                                particle_density_weight_m3
                        end
                    end
                    if !isnothing(radial_flux) ||
                       !isnothing(ionization_flux_sigma)
                        vr = radial_velocity_at_surface(
                            x, y, z, xnew, ynew, znew, vx, vy, vz,
                            flux_altitude_km[ia],
                        )
                        crossing_flux =
                            particle_density_weight_m3 * abs(vr)
                        if !isnothing(radial_flux)
                            radial_flux[ia, charge + 1, 1] += crossing_flux
                        end
                        if !isnothing(ionization_flux_sigma)
                            ionization_flux_sigma[ia, charge + 1, 1] +=
                                crossing_flux * sigma_ion_segment
                        end
                    end
                    if !isnothing(lya_rate_coefficient)
                        particle_rate_coefficient =
                            particle_density_weight_m3 * speed
                        for target in 1:3
                            lya_rate_coefficient[
                                ia, charge + 1, target, 1
                            ] += particle_rate_coefficient *
                                 sigma_lya_segment[target]
                        end
                    end
                end
            elseif altitude_after > alt
                first_surface = searchsortedlast(flux_altitude_km, alt) + 1
                last_surface = searchsortedlast(flux_altitude_km, altitude_after)
                for ia in first_surface:last_surface
                    if valid_energy
                        if !isnothing(directional_flux)
                            directional_flux[ia, ie, charge + 1, 2] +=
                                injection_flux_weight
                        end
                        if !isnothing(density_crossings)
                            density_crossings[ia, ie, charge + 1] +=
                                particle_density_weight_m3
                        end
                    end
                    if !isnothing(radial_flux) ||
                       !isnothing(ionization_flux_sigma)
                        vr = radial_velocity_at_surface(
                            x, y, z, xnew, ynew, znew, vx, vy, vz,
                            flux_altitude_km[ia],
                        )
                        crossing_flux =
                            particle_density_weight_m3 * abs(vr)
                        if !isnothing(radial_flux)
                            radial_flux[ia, charge + 1, 2] += crossing_flux
                        end
                        if !isnothing(ionization_flux_sigma)
                            ionization_flux_sigma[ia, charge + 1, 2] +=
                                crossing_flux * sigma_ion_segment
                        end
                    end
                    if !isnothing(lya_rate_coefficient)
                        particle_rate_coefficient =
                            particle_density_weight_m3 * speed
                        for target in 1:3
                            lya_rate_coefficient[
                                ia, charge + 1, target, 2
                            ] += particle_rate_coefficient *
                                 sigma_lya_segment[target]
                        end
                    end
                end
            end
            if !isnothing(surface_diagnostics) &&
               min(alt, altitude_after) <= surface_altitude_km <=
               max(alt, altitude_after)
                xc, yc, zc, _ = state_at_surface_crossing(
                    x, y, z, xnew, ynew, znew, vx, vy, vz,
                    surface_altitude_km,
                )
                position = cartesian_to_lon_lat_alt(xc, yc, zc)
                longitude = mod(position.lon_deg + 180.0, 360.0) - 180.0
                ilon = searchsortedlast(surface_longitude_edges_deg, longitude)
                ilat = searchsortedlast(
                    surface_latitude_edges_deg, position.lat_deg,
                )
                if 1 <= ilon < length(surface_longitude_edges_deg) &&
                   1 <= ilat < length(surface_latitude_edges_deg)
                    density = neutral_density(
                        model, position.lon_deg, position.lat_deg,
                        surface_altitude_km;
                        include_hot_o=cfg.include_hot_o,
                    )
                    target_density = (density.CO2, density.O, density.N2)
                    particle_rate = particle_density_weight_m3 * speed
                    # Metrics: O ionization, CO2 ionization, Ly-alpha VER,
                    # and total projectile collision-energy transfer.
                    surface_diagnostics[ilon, ilat, 1] += particle_rate *
                        target_density[2] *
                        sigma_at(model.cross_sections, charge, 2, 2, e)
                    surface_diagnostics[ilon, ilat, 2] += particle_rate *
                        target_density[1] *
                        sigma_at(model.cross_sections, charge, 1, 2, e)
                    lya_rate = 0.0
                    energy_rate = 0.0
                    for target in 1:3
                        nt = target_density[target]
                        lya_rate += nt * sigma_at(
                            model.cross_sections, charge, target, 3, e,
                        )
                        for reaction in 1:3
                            energy_rate += nt * sigma_at(
                                model.cross_sections, charge, target,
                                reaction, e,
                            ) * model.cross_sections.loss[
                                charge + 1, target, reaction
                            ]
                        end
                        energy_rate += nt * sigma_at(
                            model.cross_sections, charge, target, 4, e,
                        ) * e * elastic_loss_fraction[charge + 1, target]
                    end
                    surface_diagnostics[ilon, ilat, 3] +=
                        particle_rate * lya_rate
                    surface_diagnostics[ilon, ilat, 4] +=
                        particle_rate * energy_rate
                end
            end
        end
        x, y, z = xnew, ynew, znew
        elapsed_time += ds / speed
        tau += alpha * ds
        steps += 1
        if record
            alt_step = sqrt(x*x+y*y+z*z)/1000 - MARS_RADIUS_KM
            push!(history, HistoryEvent(id, 1, steps, collisions, elapsed_time,
                x,y,z,vx,vy,vz,alt_step,e,e,vx,vy,vz,Int8(charge),0,0,0.0))
        end
        if alpha > 0 && tau >= threshold - 8eps(threshold)
            target, reaction = choose_event(rng, model, charge, e, n)
            counts[reaction] += 1
            vx_before, vy_before, vz_before = vx, vy, vz
            theta = deg2rad(interp1(rand(rng), model.cross_sections.scatter_r,
                                    model.cross_sections.scatter_theta))
            phi = 2pi * rand(rng)
            energy_before = e
            if reaction == 4
                # Two-body elastic energy transfer in the laboratory frame.
                ratio = charge == 1 ? HP_MASS / TARGET_MASS[target] : H_MASS / TARGET_MASS[target]
                disc = max(1 - (ratio*sin(theta))^2, 0.0)
                speed_after = speed * max((ratio*cos(theta) + sqrt(disc)) / (1 + ratio), 0.0)
            else
                # Inelastic channels subtract a fixed tabulated loss.
                # State change toggles H <-> H+; Ly-alpha leaves neutral H.
                e2 = max(e - model.cross_sections.loss[charge + 1, target, reaction], 0.0)
                newcharge = reaction == 1 ? 1 - charge : (reaction == 3 ? 0 : charge)
                speed_after = sqrt(2 * e2 * QE / (newcharge == 1 ? HP_MASS : H_MASS))
                charge = newcharge
            end
            vx, vy, vz = rotate_velocity(vx, vy, vz, speed_after, theta, phi)
            collisions += 1
            if !isnothing(reaction_counts)
                collision_altitude = sqrt(x*x+y*y+z*z)/1000 - MARS_RADIUS_KM
                ibin = searchsortedlast(altitude_edges_km, collision_altitude)
                if 1 <= ibin < length(altitude_edges_km)
                    reaction_counts[ibin, reaction] += 1
                end
            end
            if record
                alt_collision = sqrt(x*x+y*y+z*z)/1000 - MARS_RADIUS_KM
                energy_after = energy(vx,vy,vz,charge)
                push!(history, HistoryEvent(id, 2, steps, collisions, elapsed_time,
                    x,y,z,vx,vy,vz,alt_collision,energy_before,energy_after,
                    vx_before,vy_before,vz_before,Int8(charge),
                    UInt8(target),UInt8(reaction),energy_before-energy_after))
            end
            # A collision ends the current free path. Sample the next one.
            tau, threshold = 0.0, -log(rand(rng))
        end
    end
    collisions == cfg.max_collisions && (stop = 5)
    final_alt = sqrt(x*x+y*y+z*z)/1000 - MARS_RADIUS_KM
    if !isnothing(thermalization_energy_map) && stop == 1 &&
       thermalization_shell_edges_km[1] <= final_alt <
       thermalization_shell_edges_km[2]
        position = cartesian_to_lon_lat_alt(x, y, z)
        longitude = mod(position.lon_deg + 180.0, 360.0) - 180.0
        ilon = searchsortedlast(surface_longitude_edges_deg, longitude)
        ilat = searchsortedlast(surface_latitude_edges_deg, position.lat_deg)
        if 1 <= ilon < length(surface_longitude_edges_deg) &&
           1 <= ilat < length(surface_latitude_edges_deg)
            shell_thickness_m = 1000.0 * (
                thermalization_shell_edges_km[2] -
                thermalization_shell_edges_km[1]
            )
            thermalization_energy_map[ilon, ilat] +=
                injection_flux_weight * energy(vx, vy, vz, charge) /
                shell_thickness_m
        end
    end
    summary = ParticleSummary(energy(vx,vy,vz,charge), final_alt, steps, collisions,
                              counts[4], counts[2], counts[3], counts[1], stop)
    if record
        push!(history, HistoryEvent(id, 3, steps, collisions, elapsed_time,
            x,y,z,vx,vy,vz,final_alt,summary.final_energy_ev,summary.final_energy_ev,
            vx,vy,vz,Int8(charge),0,0,0.0))
    end
    summary, history
end

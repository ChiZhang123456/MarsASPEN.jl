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
Transport one particle until a physical or numerical stopping condition.

Optional accumulators support three output modes without duplicating the
physics kernel:

* `record=true` stores full step and collision history.
* `reaction_counts` accumulates collision events by altitude and reaction.
* `path_length_m` accumulates `particle_weight * ds` by altitude, energy,
  and charge state.

All particles currently begin as neutral H ENA at 600 km, moving radially
toward Mars. The random stream is derived from `(cfg.seed, id)`, which makes
results reproducible and independent of thread scheduling.
"""
function run_particle_core(
    model::AspenModel, cfg::MonteCarloConfig, id::Int, record::Bool;
    altitude_edges_km::Union{Nothing,Vector{Float64}}=nothing,
    reaction_counts::Union{Nothing,Matrix{Int64}}=nothing,
    energy_edges_ev::Union{Nothing,Vector{Float64}}=nothing,
    path_length_m::Union{Nothing,Array{Float64,3}}=nothing,
)
    # A separate deterministic RNG per particle prevents thread-order effects.
    rng = Xoshiro(hash((cfg.seed, id)))
    x, y, z = (MARS_RADIUS_KM + cfg.initial_altitude_km) * 1000, 0.0, 0.0
    vx, vy, vz = -cfg.initial_speed_m_s, 0.0, 0.0
    # charge=0 is neutral H ENA. `tau` is accumulated optical depth and
    # `threshold` is the exponentially distributed optical depth to collision.
    charge, tau, threshold = 0, 0.0, -log(rand(rng))
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
                path_length_m[ia, ie, charge + 1] += cfg.particle_weight * ds
            end
        end
        x += vx / speed * ds; y += vy / speed * ds; z += vz / speed * ds
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
    summary = ParticleSummary(energy(vx,vy,vz,charge), final_alt, steps, collisions,
                              counts[4], counts[2], counts[3], counts[1], stop)
    if record
        push!(history, HistoryEvent(id, 3, steps, collisions, elapsed_time,
            x,y,z,vx,vy,vz,final_alt,summary.final_energy_ev,summary.final_energy_ev,
            vx,vy,vz,Int8(charge),0,0,0.0))
    end
    summary, history
end

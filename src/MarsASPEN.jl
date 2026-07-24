module MarsASPEN

using MAT
using Random
using Statistics
using Base.Threads

export AspenModel, MonteCarloConfig, ParticleSummary, HistoryEvent, load_model,
       run_ensemble, run_detailed_ensemble, write_detailed_mat

const QE = 1.602176634e-19
const AMU = 1.66053906660e-27
const MARS_RADIUS_KM = 3388.25
const H_MASS = 1.00782503223 * AMU
const HP_MASS = 1.007276466621 * AMU
const TARGET_MASS = (44.0095 * AMU, 15.999 * AMU, 28.0134 * AMU)
const ATMOSPHERE_MASS = (44.01 * AMU, 15.999 * AMU, 28.014 * AMU)
const KB = 1.380649e-23
const MARS_G0 = 3.71
const REACTION_NAMES = (:state_change, :ionization, :lya, :elastic)

struct Atmosphere
    lon::Vector{Float64}
    lat::Vector{Float64}
    alt::Vector{Float64}
    logn::Array{Float64,4} # lon, lat, alt, target
    tn::Array{Float64,3}
end

struct CrossSections
    energy::Vector{Float64}
    sigma::Array{Float64,4} # charge+1, target, reaction, energy
    loss::Array{Float64,3}
    scatter_r::Vector{Float64}
    scatter_theta::Vector{Float64}
end

struct HotAtmosphere
    lon::Vector{Float64}
    lat::Vector{Float64}
    alt::Vector{Float64}
    logn_o::Array{Float64,3}
end

struct AspenModel
    atmosphere::Atmosphere
    hot_atmosphere::HotAtmosphere
    cross_sections::CrossSections
end

Base.@kwdef struct MonteCarloConfig
    n_particles::Int = 1000
    initial_altitude_km::Float64 = 600.0
    initial_speed_m_s::Float64 = 400_000.0
    seed::Int = 7
    safety_factor::Float64 = 0.4
    max_step_m::Float64 = 1000.0
    min_energy_ev::Float64 = 10.0
    min_altitude_km::Float64 = 100.0
    max_altitude_km::Float64 = 1000.0
    max_collisions::Int = 2000
    max_steps_per_collision::Int = 100_000
    include_hot_o::Bool = true
end

struct ParticleSummary
    final_energy_ev::Float64
    final_altitude_km::Float64
    n_steps::Int
    n_collisions::Int
    n_elastic::Int
    n_ionization::Int
    n_lya::Int
    n_state_change::Int
    stop_code::UInt8
end

struct HistoryEvent
    particle_id::Int
    event_type::UInt8 # 0 initial, 1 transport, 2 collision, 3 final
    step::Int
    collision::Int
    time_s::Float64
    x_m::Float64
    y_m::Float64
    z_m::Float64
    vx_m_s::Float64
    vy_m_s::Float64
    vz_m_s::Float64
    altitude_km::Float64
    energy_before_ev::Float64
    energy_ev::Float64
    vx_before_m_s::Float64
    vy_before_m_s::Float64
    vz_before_m_s::Float64
    charge_state::Int8
    target::UInt8 # 0 none, 1 CO2, 2 O, 3 N2
    reaction::UInt8 # 0 none, 1 state change, 2 ionization, 3 Ly-alpha, 4 elastic
    energy_loss_ev::Float64
end

function _read_table(path::AbstractString, nskip::Int)
    rows = Vector{Vector{Float64}}()
    for (i, line) in enumerate(eachline(path))
        i <= nskip && continue
        s = strip(line)
        (isempty(s) || startswith(s, "#")) && continue
        push!(rows, parse.(Float64, split(s)))
    end
    reduce(vcat, permutedims.(rows))
end

@inline function interp1(x::Float64, xp::Vector{Float64}, fp)
    (x < xp[1] || x > xp[end]) && return 0.0
    i = clamp(searchsortedlast(xp, x), 1, length(xp) - 1)
    w = (x - xp[i]) / (xp[i + 1] - xp[i])
    muladd(w, fp[i + 1] - fp[i], fp[i])
end

function _parse_losses(path::AbstractString, charge::Int)
    names = charge == 1 ?
        Dict("sigma_10"=>1, "sigma_ip"=>2, "sigma_La"=>3, "sigma_el"=>4) :
        Dict("sigma_01"=>1, "sigma_ia"=>2, "sigma_La"=>3, "sigma_el"=>4)
    out = zeros(4)
    for line in eachline(path)
        !startswith(line, "#") && break
        for (name, ir) in names
            occursin(name, line) || continue
            m = match(r"Q\s*=\s*([+-]?\d+(?:\.\d+)?)\s*eV", line)
            m === nothing || (out[ir] = abs(parse(Float64, m.captures[1])))
        end
    end
    out
end

function load_model(repo_root::AbstractString=pkgdir(@__MODULE__);
                    ls::Int=0, f107::Int=70,
                    atmosphere_data_dir::Union{Nothing,AbstractString}=nothing)
    atmosphere_dir = isnothing(atmosphere_data_dir) ?
        get(ENV, "MARSASPEN_ATMOSPHERE_DIR", joinpath(repo_root, "data", "atmosphere")) :
        String(atmosphere_data_dir)
    matpath = joinpath(atmosphere_dir,
                       "gitm_ls$(lpad(ls, 3, '0'))_f$(lpad(f107, 3, '0')).mat")
    isfile(matpath) || error(
        "MGITM file not found: $matpath. Set MARSASPEN_ATMOSPHERE_DIR or " *
        "pass atmosphere_data_dir to load_model.",
    )
    d = matread(matpath)
    lon, lat, alt = vec(d["lon_deg"]), vec(d["lat_deg"]), vec(d["alt_km"])
    logn = Array{Float64}(undef, length(lon), length(lat), length(alt), 3)
    for (it, field) in enumerate(("nCO2", "nO", "nN2"))
        logn[:,:,:,it] .= log.(max.(Float64.(d[field]), 1e-300))
    end
    atmosphere = Atmosphere(lon, lat, alt, logn, Float64.(d["Tn"]))
    hotpath = joinpath(atmosphere_dir,
                       "mamps_ls$(lpad(ls, 3, '0'))_f$(lpad(f107, 3, '0')).mat")
    isfile(hotpath) || error("MAMPS file not found: $hotpath")
    hot = matread(hotpath)
    hot_atmosphere = HotAtmosphere(
        vec(hot["lon_deg"]), vec(hot["lat_deg"]), vec(hot["alt_km"]),
        log.(max.(Float64.(hot["nO_hot"]), 1e-300)),
    )

    csdir = joinpath(repo_root, "data", "cross_sections")
    files = (("Hplus_CO2_cross_sections.txt", 1, 1),
             ("Hplus_O_cross_sections.txt", 1, 2),
             ("Hplus_N2_cross_sections.txt", 1, 3),
             ("H_CO2_cross_sections.txt", 0, 1),
             ("H_O_cross_sections.txt", 0, 2),
             ("H_N2_cross_sections.txt", 0, 3))
    energy = 10.0 .^ range(0.0, 7.0, length=768)
    sigma = zeros(2, 3, 4, length(energy))
    loss = zeros(2, 3, 4)
    for (name, charge, target) in files
        path = joinpath(csdir, name)
        raw = _read_table(path, 9)
        columns = charge == 1 ? (2, 3, 4, 5) : (3, 5, 4, 2)
        for ir in 1:4, ie in eachindex(energy)
            sigma[charge + 1, target, ir, ie] =
                interp1(energy[ie], vec(raw[:,1]), view(raw,:,columns[ir])) * 1e-4
        end
        loss[charge + 1, target, :] .= _parse_losses(path, charge)
    end
    scatter = _read_table(joinpath(csdir, "scattering_angle_distribution.txt"), 7)
    cross_sections = CrossSections(energy, sigma, loss, vec(scatter[:,1]), vec(scatter[:,2]))
    AspenModel(atmosphere, hot_atmosphere, cross_sections)
end

@inline function bracket(grid::Vector{Float64}, x::Float64)
    i = clamp(searchsortedlast(grid, x), 1, length(grid) - 1)
    w = clamp((x - grid[i]) / (grid[i + 1] - grid[i]), 0.0, 1.0)
    i, i + 1, w
end

@inline function hot_o_density(a::HotAtmosphere, lon::Float64, lat::Float64, alt::Float64)
    (alt < a.alt[1] || alt > a.alt[end]) && return 0.0
    i0, i1, wx = lonbracket(a.lon, lon)
    j0, j1, wy = bracket(a.lat, lat)
    k0, k1, wz = bracket(a.alt, alt)
    c00 = muladd(wx, a.logn_o[i1,j0,k0] - a.logn_o[i0,j0,k0], a.logn_o[i0,j0,k0])
    c10 = muladd(wx, a.logn_o[i1,j1,k0] - a.logn_o[i0,j1,k0], a.logn_o[i0,j1,k0])
    c01 = muladd(wx, a.logn_o[i1,j0,k1] - a.logn_o[i0,j0,k1], a.logn_o[i0,j0,k1])
    c11 = muladd(wx, a.logn_o[i1,j1,k1] - a.logn_o[i0,j1,k1], a.logn_o[i0,j1,k1])
    c0 = muladd(wy, c10 - c00, c00)
    c1 = muladd(wy, c11 - c01, c01)
    exp(muladd(wz, c1 - c0, c0))
end

@inline function lonbracket(grid::Vector{Float64}, xraw::Float64)
    x = mod(xraw - grid[1], 360.0) + grid[1]
    i = searchsortedlast(grid, x)
    if i == 0
        return length(grid), 1, (x + 360.0 - grid[end]) / (grid[1] + 360.0 - grid[end])
    elseif i == length(grid)
        return i, 1, (x - grid[i]) / (grid[1] + 360.0 - grid[i])
    end
    i, i + 1, (x - grid[i]) / (grid[i + 1] - grid[i])
end

@inline function density3(a::Atmosphere, lon::Float64, lat::Float64, alt::Float64)
    i0, i1, wx = lonbracket(a.lon, lon)
    j0, j1, wy = bracket(a.lat, lat)
    k0, k1, wz = bracket(a.alt, clamp(alt, a.alt[1], a.alt[end]))
    density = ntuple(3) do it
        c00 = muladd(wx, a.logn[i1,j0,k0,it] - a.logn[i0,j0,k0,it], a.logn[i0,j0,k0,it])
        c10 = muladd(wx, a.logn[i1,j1,k0,it] - a.logn[i0,j1,k0,it], a.logn[i0,j1,k0,it])
        c01 = muladd(wx, a.logn[i1,j0,k1,it] - a.logn[i0,j0,k1,it], a.logn[i0,j0,k1,it])
        c11 = muladd(wx, a.logn[i1,j1,k1,it] - a.logn[i0,j1,k1,it], a.logn[i0,j1,k1,it])
        c0 = muladd(wy, c10 - c00, c00)
        c1 = muladd(wy, c11 - c01, c01)
        exp(muladd(wz, c1 - c0, c0))
    end
    if alt > a.alt[end]
        t00 = muladd(wx, a.tn[i1,j0,k1] - a.tn[i0,j0,k1], a.tn[i0,j0,k1])
        t10 = muladd(wx, a.tn[i1,j1,k1] - a.tn[i0,j1,k1], a.tn[i0,j1,k1])
        tn_top = muladd(wy, t10 - t00, t00)
        g = MARS_G0 * (MARS_RADIUS_KM / (MARS_RADIUS_KM + a.alt[end]))^2
        return ntuple(it -> density[it] * exp(-(alt - a.alt[end]) /
            (KB * tn_top / (ATMOSPHERE_MASS[it] * g) / 1000)), 3)
    end
    density
end

@inline energy(vx, vy, vz, charge) =
    0.5 * (charge == 1 ? HP_MASS : H_MASS) * (vx*vx + vy*vy + vz*vz) / QE

@inline function sigma_at(cs::CrossSections, charge, target, reaction, e)
    interp1(e, cs.energy, view(cs.sigma, charge + 1, target, reaction, :))
end

@inline function local_state(model::AspenModel, x, y, z, charge, e, include_hot_o::Bool)
    rkm = sqrt(x*x + y*y + z*z) / 1000
    alt = rkm - MARS_RADIUS_KM
    lon = mod(rad2deg(atan(y, x)), 360.0)
    lat = rad2deg(asin(clamp(z / (rkm * 1000), -1.0, 1.0)))
    n = density3(model.atmosphere, lon, lat, alt)
    if include_hot_o
        n = (n[1], n[2] + hot_o_density(model.hot_atmosphere, lon, lat, alt), n[3])
    end
    alpha = 0.0
    @inbounds for it in 1:3, ir in 1:4
        alpha += n[it] * sigma_at(model.cross_sections, charge, it, ir, e)
    end
    alt, n, alpha
end

@inline function choose_event(rng, model, charge, e, n)
    weights = ntuple(12) do q
        it = (q - 1) ÷ 4 + 1
        ir = (q - 1) % 4 + 1
        n[it] * sigma_at(model.cross_sections, charge, it, ir, e)
    end
    total = sum(weights)
    u = rand(rng) * total
    acc = 0.0
    @inbounds for q in 1:12
        acc += weights[q]
        if u <= acc
            return (q - 1) ÷ 4 + 1, (q - 1) % 4 + 1
        end
    end
    3, 4
end

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

function run_particle_core(model::AspenModel, cfg::MonteCarloConfig, id::Int, record::Bool)
    rng = Xoshiro(hash((cfg.seed, id)))
    x, y, z = (MARS_RADIUS_KM + cfg.initial_altitude_km) * 1000, 0.0, 0.0
    vx, vy, vz = -cfg.initial_speed_m_s, 0.0, 0.0
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
        candidate = alpha > 0 ? min(cfg.safety_factor / alpha, cfg.max_step_m) : cfg.max_step_m
        remaining = threshold - tau
        ds = alpha > 0 ? min(candidate, remaining / alpha) : candidate
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
                ratio = charge == 1 ? HP_MASS / TARGET_MASS[target] : H_MASS / TARGET_MASS[target]
                disc = max(1 - (ratio*sin(theta))^2, 0.0)
                speed_after = speed * max((ratio*cos(theta) + sqrt(disc)) / (1 + ratio), 0.0)
            else
                e2 = max(e - model.cross_sections.loss[charge + 1, target, reaction], 0.0)
                newcharge = reaction == 1 ? 1 - charge : (reaction == 3 ? 0 : charge)
                speed_after = sqrt(2 * e2 * QE / (newcharge == 1 ? HP_MASS : H_MASS))
                charge = newcharge
            end
            vx, vy, vz = rotate_velocity(vx, vy, vz, speed_after, theta, phi)
            collisions += 1
            if record
                alt_collision = sqrt(x*x+y*y+z*z)/1000 - MARS_RADIUS_KM
                energy_after = energy(vx,vy,vz,charge)
                push!(history, HistoryEvent(id, 2, steps, collisions, elapsed_time,
                    x,y,z,vx,vy,vz,alt_collision,energy_before,energy_after,
                    vx_before,vy_before,vz_before,Int8(charge),
                    UInt8(target),UInt8(reaction),energy_before-energy_after))
            end
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

function run_detailed_ensemble(model::AspenModel, cfg::MonteCarloConfig)
    summaries = Vector{ParticleSummary}(undef, cfg.n_particles)
    histories = Vector{Vector{HistoryEvent}}(undef, cfg.n_particles)
    @threads for i in 1:cfg.n_particles
        summaries[i], histories[i] = run_particle_core(model, cfg, i, true)
    end
    summaries, histories
end

function write_detailed_mat(filename::AbstractString, summaries, histories;
                            config::MonteCarloConfig)
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
        "config_max_step_m" => config.max_step_m,
    )
    mkpath(dirname(abspath(filename)))
    matwrite(filename, data; compress=true)
    filename
end

end

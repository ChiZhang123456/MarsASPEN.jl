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

function load_cross_sections(repo_root::AbstractString)
    cross_section_dir = joinpath(repo_root, "data", "cross_sections")
    files = (
        ("Hplus_CO2_cross_sections.txt", 1, 1),
        ("Hplus_O_cross_sections.txt", 1, 2),
        ("Hplus_N2_cross_sections.txt", 1, 3),
        ("H_CO2_cross_sections.txt", 0, 1),
        ("H_O_cross_sections.txt", 0, 2),
        ("H_N2_cross_sections.txt", 0, 3),
    )
    energy = 10.0 .^ range(0.0, 7.0, length=768)
    sigma = zeros(2, 3, 4, length(energy))
    loss = zeros(2, 3, 4)
    for (name, charge, target) in files
        path = joinpath(cross_section_dir, name)
        raw = _read_table(path, 9)
        columns = charge == 1 ? (2, 3, 4, 5) : (3, 5, 4, 2)
        for reaction in 1:4, ie in eachindex(energy)
            sigma[charge + 1, target, reaction, ie] =
                interp1(
                    energy[ie], vec(raw[:,1]), view(raw,:,columns[reaction]),
                ) * 1e-4
        end
        loss[charge + 1, target, :] .= _parse_losses(path, charge)
    end
    scatter = _read_table(
        joinpath(cross_section_dir, "scattering_angle_distribution.txt"), 7,
    )
    CrossSections(
        energy, sigma, loss, vec(scatter[:,1]), vec(scatter[:,2]),
    )
end

@inline function sigma_at(cs::CrossSections, charge, target, reaction, e)
    interp1(e, cs.energy, view(cs.sigma, charge + 1, target, reaction, :))
end

@inline function local_state(model::AspenModel, x, y, z, charge, e, include_hot_o::Bool)
    rkm = sqrt(x*x + y*y + z*z) / 1000
    alt = rkm - MARS_RADIUS_KM
    lon = mod(rad2deg(atan(y, x)), 360.0)
    lat = rad2deg(asin(clamp(z / (rkm * 1000), -1.0, 1.0)))
    n = transport_density3(model.atmosphere, lon, lat, alt)
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

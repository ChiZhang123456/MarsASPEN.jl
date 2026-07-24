function normalize_f107(solar)
    key = lowercase(replace(string(solar), r"[_\-\s]" => ""))
    key in ("solarmax", "max", "f200", "200") && return 200
    key in ("solarmoderate", "moderate", "mod", "solarmean", "mean", "f130", "130") && return 130
    key in ("solarmin", "min", "f070", "f70", "70") && return 70
    throw(ArgumentError("solar must be solar_max, solar_moderate, solar_min, 200, 130, or 70"))
end

function available_atmosphere_cases()
    [(ls=ls, f107=f107) for ls in (0, 90, 180, 270), f107 in (70, 130, 200)] |> vec
end

function load_atmospheres(
    atmosphere_dir::AbstractString, ls::Int, f107_value::Int,
)
    matpath = joinpath(atmosphere_dir,
                       "gitm_ls$(lpad(ls, 3, '0'))_f$(lpad(f107_value, 3, '0')).mat")
    isfile(matpath) || error(
        "MGITM file not found: $matpath. Set MARSASPEN_ATMOSPHERE_DIR or " *
        "pass atmosphere_data_dir to load_model.",
    )
    d = matread(matpath)
    lon, lat, alt = vec(d["lon_deg"]), vec(d["lat_deg"]), vec(d["alt_km"])
    logn = Array{Float64}(undef, length(lon), length(lat), length(alt), 5)
    for (it, field) in enumerate(("nCO2", "nO", "nO2", "nN2", "nCO"))
        logn[:,:,:,it] .= log.(max.(Float64.(d[field]), 1e-300))
    end
    atmosphere = Atmosphere(lon, lat, alt, logn, Float64.(d["Tn"]))
    hotpath = joinpath(atmosphere_dir,
                       "mamps_ls$(lpad(ls, 3, '0'))_f$(lpad(f107_value, 3, '0')).mat")
    isfile(hotpath) || error("MAMPS file not found: $hotpath")
    hot = matread(hotpath)
    hot_atmosphere = HotAtmosphere(
        vec(hot["lon_deg"]), vec(hot["lat_deg"]), vec(hot["alt_km"]),
        log.(max.(Float64.(hot["nO_hot"]), 1e-300)),
    )
    atmosphere, hot_atmosphere
end

function load_model(repo_root::AbstractString=pkgdir(@__MODULE__);
                    ls::Int=0, solar=70, f107::Union{Nothing,Int}=nothing,
                    atmosphere_data_dir::Union{Nothing,AbstractString}=nothing)
    ls in (0, 90, 180, 270) ||
        throw(ArgumentError("ls must be one of 0, 90, 180, or 270"))
    f107_value = isnothing(f107) ? normalize_f107(solar) : normalize_f107(f107)
    atmosphere_dir = isnothing(atmosphere_data_dir) ?
        get(ENV, "MARSASPEN_ATMOSPHERE_DIR", joinpath(repo_root, "data", "atmosphere")) :
        String(atmosphere_data_dir)
    atmosphere, hot_atmosphere =
        load_atmospheres(atmosphere_dir, ls, f107_value)
    cross_sections = load_cross_sections(repo_root)
    AspenModel(atmosphere, hot_atmosphere, cross_sections)
end

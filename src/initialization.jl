# Model initialization. The original processed GITM and AMPS MATLAB files are
# read directly. Their positions are already expressed in MSO coordinates.

"""Normalize a solar-activity label to the supported F10.7 value."""
function normalize_f107(solar)
    key = lowercase(replace(string(solar), r"[_\-\s]" => ""))
    key in ("solarmax", "max", "f200", "200") && return 200
    key in ("solarmoderate", "moderate", "mod", "solarmean", "mean", "f130", "130") && return 130
    key in ("solarmin", "min", "f070", "f70", "70") && return 70
    throw(ArgumentError("solar must be solar_max, solar_moderate, solar_min, 200, 130, or 70"))
end

"""Return all available `(Ls, F10.7)` neutral-atmosphere combinations."""
function available_atmosphere_cases()
    [(ls=ls, f107=f107) for ls in (0, 90, 180, 270), f107 in (70, 130, 200)] |> vec
end

"""Return the original GITM and AMPS case names for one seasonal condition."""
function raw_atmosphere_case(ls::Int)
    ls == 0   && return ("aequ", "subL0")
    ls == 90  && return ("aph",  "subL0")
    ls == 180 && return ("aequ", "subL180")
    ls == 270 && return ("per",  "subL270")
    throw(ArgumentError("ls must be one of 0, 90, 180, or 270"))
end

"""Read one original vector-form GITM MAT file without coordinate rotation."""
function read_raw_gitm(path::AbstractString)
    d = matread(path)
    nlon_source, nlat, nalt = 73, 36, 50
    cube(field) = reshape(
        vec(Float64.(d[field])), nlon_source, nlat, nalt,
    )[1:72, :, :]
    lon = vec(cube("Longitude")[:, 1, 1])
    lat = vec(cube("Latitude")[1, :, 1])
    alt = vec(cube("Altitude")[1, 1, :])
    (
        lon=lon,
        lat=lat,
        alt=alt,
        temperature=cube("temp"),
        # Original densities are cm^-3. MarsASPEN uses m^-3.
        nco2=cube("NCO2") .* 1e6,
        no=cube("NO") .* 1e6,
    )
end

"""Read one original vector-form AMPS MAT file without coordinate rotation."""
function read_raw_amps(path::AbstractString)
    d = matread(path)
    nlon_source, nlat, nalt = 73, 36, 167
    cube(field) = reshape(
        vec(Float64.(d[field])), nlon_source, nlat, nalt,
    )[1:72, :, :]
    altitude_cube = cube("Altitude")
    (
        lon=vec(cube("Longitude")[:, 1, 1]),
        lat=vec(cube("Latitude")[1, :, 1]),
        # Geometric altitude differs by less than a metre within a shell.
        alt=vec(dropdims(mean(altitude_cube; dims=(1, 2)); dims=(1, 2))),
        # Original hot-O density is cm^-3.
        no_hot=cube("dens_oh") .* 1e6,
    )
end

"""Interpolate positive density fields logarithmically."""
@inline log_blend(low, high, weight) =
    exp.((1 - weight) .* log.(max.(low, 1e-300)) .+
         weight .* log.(max.(high, 1e-300)))

"""
Load original processed GITM and AMPS files from `GITM/` and `AMPS/`.

Solar minimum and maximum are read directly. The source collection has no
separate moderate case, so F10.7 = 130 is evaluated between F10.7 = 70 and
200. Density is interpolated in log space and temperature linearly.
"""
function load_atmospheres(
    atmosphere_dir::AbstractString, ls::Int, f107_value::Int,
)
    season, subl = raw_atmosphere_case(ls)
    function paths(activity)
        (
            joinpath(
                atmosphere_dir, "GITM",
                "gitm_$(season)$(activity)_$(subl)_alt220.mat",
            ),
            joinpath(
                atmosphere_dir, "AMPS",
                "dsmc_$(season)$(activity).mat",
            ),
        )
    end
    low_paths = paths("min")
    high_paths = paths("max")
    for path in (low_paths..., high_paths...)
        isfile(path) || error(
            "Atmosphere file not found: $path. Set MARSASPEN_ATMOSPHERE_DIR " *
            "or pass atmosphere_data_dir to load_model.",
        )
    end
    if f107_value == 70
        cold = read_raw_gitm(low_paths[1])
        hot = read_raw_amps(low_paths[2])
    elseif f107_value == 200
        cold = read_raw_gitm(high_paths[1])
        hot = read_raw_amps(high_paths[2])
    else
        weight = (f107_value - 70) / (200 - 70)
        cold_low, cold_high =
            read_raw_gitm(low_paths[1]), read_raw_gitm(high_paths[1])
        hot_low, hot_high =
            read_raw_amps(low_paths[2]), read_raw_amps(high_paths[2])
        cold = (
            lon=cold_low.lon, lat=cold_low.lat, alt=cold_low.alt,
            temperature=(1 - weight) .* cold_low.temperature .+
                        weight .* cold_high.temperature,
            nco2=log_blend(cold_low.nco2, cold_high.nco2, weight),
            no=log_blend(cold_low.no, cold_high.no, weight),
        )
        hot = (
            lon=hot_low.lon, lat=hot_low.lat, alt=hot_low.alt,
            no_hot=log_blend(hot_low.no_hot, hot_high.no_hot, weight),
        )
    end
    logn = fill(log(1e-300), length(cold.lon), length(cold.lat),
                length(cold.alt), 5)
    logn[:, :, :, 1] .= log.(max.(cold.nco2, 1e-300))
    logn[:, :, :, 2] .= log.(max.(cold.no, 1e-300))
    atmosphere = Atmosphere(
        cold.lon, cold.lat, cold.alt, logn, cold.temperature,
    )
    hot_atmosphere = HotAtmosphere(
        hot.lon, hot.lat, hot.alt, log.(max.(hot.no_hot, 1e-300)),
    )
    atmosphere, hot_atmosphere
end

"""
Construct a complete transport model.

Keyword arguments select season and solar activity. `atmosphere_data_dir`
allows the large atmosphere files to live outside the package checkout.
"""
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

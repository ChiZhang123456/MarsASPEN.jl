# Neutral-atmosphere interpolation and coordinate conversion.
#
# Horizontal interpolation is periodic in longitude and linear in latitude.
# Density interpolation is linear in log density. From 98.75 down to 80 km the
# lowest two MGITM layers are log-linearly extrapolated. Above the MGITM top,
# each cold species follows a hydrostatic exponential using the local top-layer
# temperature. MAMPS hot O is used only inside its native altitude range.

"""Return bracketing grid indices and linear interpolation weight."""
@inline function bracket(grid::Vector{Float64}, x::Float64)
    i = clamp(searchsortedlast(grid, x), 1, length(grid) - 1)
    w = clamp((x - grid[i]) / (grid[i + 1] - grid[i]), 0.0, 1.0)
    i, i + 1, w
end
"""Altitude bracket that permits MGITM extrapolation down to 80 km."""
@inline function atmosphere_altitude_bracket(grid::Vector{Float64}, altitude::Float64)
    effective_altitude = max(altitude, MODEL_MIN_ALTITUDE_KM)
    if effective_altitude < grid[1]
        return 1, 2, (effective_altitude - grid[1]) / (grid[2] - grid[1])
    end
    bracket(grid, effective_altitude)
end

"""Trilinearly interpolate MAMPS hot-O number density in m^-3."""
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

"""Return periodic longitude indices and interpolation weight."""
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

"""
Interpolate all five MGITM species at one longitude, latitude, and altitude.

The returned tuple order is `(CO2, O, O2, N2, CO)`, all in m^-3.
"""
@inline function density3(a::Atmosphere, lon::Float64, lat::Float64, alt::Float64)
    i0, i1, wx = lonbracket(a.lon, lon)
    j0, j1, wy = bracket(a.lat, lat)
    k0, k1, wz = atmosphere_altitude_bracket(a.alt, min(alt, a.alt[end]))
    density = ntuple(5) do it
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
            (KB * tn_top / (ATMOSPHERE_MASS[it] * g) / 1000)), 5)
    end
    density
end

"""
Fast density path used by transport.

Only CO2, O, and N2 are returned because the present collision database does
not contain O2 or CO projectile cross sections.
"""
@inline function transport_density3(
    a::Atmosphere, lon::Float64, lat::Float64, alt::Float64,
)
    i0, i1, wx = lonbracket(a.lon, lon)
    j0, j1, wy = bracket(a.lat, lat)
    k0, k1, wz = atmosphere_altitude_bracket(a.alt, min(alt, a.alt[end]))
    field_indices = (1, 2, 4)
    density = ntuple(3) do target
        it = field_indices[target]
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
        mass_indices = (1, 2, 4)
        return ntuple(target -> density[target] * exp(-(alt - a.alt[end]) /
            (KB * tn_top / (ATMOSPHERE_MASS[mass_indices[target]] * g) / 1000)), 3)
    end
    density
end

"""Interpolate MGITM neutral temperature in K."""
@inline function temperature3(a::Atmosphere, lon::Float64, lat::Float64, alt::Float64)
    i0, i1, wx = lonbracket(a.lon, lon)
    j0, j1, wy = bracket(a.lat, lat)
    k0, k1, wz = atmosphere_altitude_bracket(a.alt, min(alt, a.alt[end]))
    c00 = muladd(wx, a.tn[i1,j0,k0] - a.tn[i0,j0,k0], a.tn[i0,j0,k0])
    c10 = muladd(wx, a.tn[i1,j1,k0] - a.tn[i0,j1,k0], a.tn[i0,j1,k0])
    c01 = muladd(wx, a.tn[i1,j0,k1] - a.tn[i0,j0,k1], a.tn[i0,j0,k1])
    c11 = muladd(wx, a.tn[i1,j1,k1] - a.tn[i0,j1,k1], a.tn[i0,j1,k1])
    c0 = muladd(wy, c10 - c00, c00)
    c1 = muladd(wy, c11 - c01, c01)
    muladd(wz, c1 - c0, c0)
end

"""Return named neutral densities and temperature at a geographic position."""
function neutral_density(model::AspenModel, lon_deg::Real, lat_deg::Real,
                         altitude_km::Real; include_hot_o::Bool=true)
    cold = density3(
        model.atmosphere, Float64(lon_deg), Float64(lat_deg), Float64(altitude_km),
    )
    hot = include_hot_o ?
        hot_o_density(
            model.hot_atmosphere, Float64(lon_deg), Float64(lat_deg), Float64(altitude_km),
        ) : 0.0
    (
        CO2=cold[1], O=cold[2] + hot, O2=cold[3], N2=cold[4], CO=cold[5],
        Tn=temperature3(
            model.atmosphere, Float64(lon_deg), Float64(lat_deg), Float64(altitude_km),
        ),
        O_cold=cold[2], O_hot=hot,
    )
end

"""Convert Mars-centered Cartesian coordinates to longitude, latitude, altitude."""
function cartesian_to_lon_lat_alt(x::Real, y::Real, z::Real;
                                  position_unit::Symbol=:m)
    scale = position_unit === :m ? 1 / 1000 :
            position_unit === :km ? 1.0 :
            throw(ArgumentError("position_unit must be :m or :km"))
    xkm, ykm, zkm = Float64(x) * scale, Float64(y) * scale, Float64(z) * scale
    radius = sqrt(xkm*xkm + ykm*ykm + zkm*zkm)
    radius > 0 || throw(ArgumentError("Cartesian position must have nonzero radius"))
    (
        lon_deg=mod(rad2deg(atan(ykm, xkm)), 360.0),
        lat_deg=rad2deg(asin(clamp(zkm / radius, -1.0, 1.0))),
        altitude_km=radius - MARS_RADIUS_KM,
    )
end

"""Evaluate `neutral_density` from a Mars-centered Cartesian position."""
function neutral_density_xyz(model::AspenModel, x::Real, y::Real, z::Real;
                             position_unit::Symbol=:m, include_hot_o::Bool=true)
    position = cartesian_to_lon_lat_alt(x, y, z; position_unit=position_unit)
    neutral_density(
        model, position.lon_deg, position.lat_deg, position.altitude_km;
        include_hot_o=include_hot_o,
    )
end

@inline energy(vx, vy, vz, charge) =
    0.5 * (charge == 1 ? HP_MASS : H_MASS) * (vx*vx + vy*vy + vz*vz) / QE

# Three-dimensional spherical-grid diagnostics and MAT serialization.
#
# The transport step is at most 1 km in the production examples. Each segment
# is therefore assigned by its midpoint to a 1 km altitude cell and to the
# matching MGITM horizontal cell. Track residence time gives number density,
# path length gives total scalar flux, and radial displacement gives signed
# and directional radial flux.

"""Create a zeroed spherical grid using the native MGITM horizontal spacing."""
function create_spatial_grid(
    model::AspenModel;
    altitude_edges_km::AbstractVector{<:Real}=collect(80.0:1.0:600.0),
)
    lon_centers = model.atmosphere.lon
    lat_centers = model.atmosphere.lat
    dlon = median(diff(lon_centers))
    dlat = median(diff(lat_centers))
    all(isapprox.(diff(lon_centers), dlon; rtol=0, atol=1e-10)) ||
        throw(ArgumentError("MGITM longitude grid must have uniform spacing"))
    all(isapprox.(diff(lat_centers), dlat; rtol=0, atol=1e-10)) ||
        throw(ArgumentError("MGITM latitude grid must have uniform spacing"))
    lon_edges = collect(
        (lon_centers[1] - 0.5dlon):dlon:
        (lon_centers[end] + 0.5dlon),
    )
    lat_edges = collect(
        (lat_centers[1] - 0.5dlat):dlat:
        (lat_centers[end] + 0.5dlat),
    )
    lat_edges[1], lat_edges[end] = -90.0, 90.0
    alt_edges = Float64.(altitude_edges_km)
    all(diff(alt_edges) .> 0) ||
        throw(ArgumentError("altitude_edges_km must increase"))
    nlon, nlat, nalt =
        length(lon_edges)-1, length(lat_edges)-1, length(alt_edges)-1
    volume = zeros(Float64, nlat, nalt)
    delta_lon_rad = deg2rad(dlon)
    for ilat in 1:nlat, ialt in 1:nalt
        solid_angle = delta_lon_rad * (
            sin(deg2rad(lat_edges[ilat + 1])) -
            sin(deg2rad(lat_edges[ilat]))
        )
        r1 = (MARS_RADIUS_KM + alt_edges[ialt]) * 1000
        r2 = (MARS_RADIUS_KM + alt_edges[ialt + 1]) * 1000
        volume[ilat, ialt] = (r2^3 - r1^3) * solid_angle / 3
    end
    dims_charge = (nlon, nlat, nalt, 2)
    SpatialGridAccumulator(
        lon_edges, lat_edges, alt_edges, volume,
        zeros(Float64, dims_charge),
        zeros(Float64, dims_charge),
        zeros(Float64, dims_charge),
        zeros(Float64, dims_charge),
        zeros(Float64, dims_charge),
        zeros(Float64, nlon, nlat, nalt, 2, 3, 4),
        zeros(UInt32, nlon, nlat, nalt, 4),
        zeros(Float64, dims_charge),
        zeros(Float64, dims_charge),
        [ReentrantLock() for _ in 1:nlon],
    )
end

"""Return spherical-grid indices and cell volume, or `nothing` if outside."""
@inline function spatial_grid_cell(
    grid::SpatialGridAccumulator, x::Float64, y::Float64, z::Float64,
)
    radius_m = sqrt(x*x + y*y + z*z)
    altitude_km = radius_m / 1000 - MARS_RADIUS_KM
    longitude_deg = mod(rad2deg(atan(y, x)), 360.0)
    latitude_deg = rad2deg(asin(clamp(z / radius_m, -1.0, 1.0)))
    ilon = searchsortedlast(grid.longitude_edges_deg, longitude_deg)
    ilat = searchsortedlast(grid.latitude_edges_deg, latitude_deg)
    ialt = searchsortedlast(grid.altitude_edges_km, altitude_km)
    if 1 <= ilon < length(grid.longitude_edges_deg) &&
       1 <= ilat < length(grid.latitude_edges_deg) &&
       1 <= ialt < length(grid.altitude_edges_km)
        return ilon, ilat, ialt, grid.cell_volume_m3[ilat, ialt]
    end
    nothing
end

"""Accumulate residence density and flux moments for one free-flight segment."""
@inline function accumulate_spatial_segment!(
    grid::SpatialGridAccumulator,
    x::Float64, y::Float64, z::Float64,
    xnew::Float64, ynew::Float64, znew::Float64,
    vx::Float64, vy::Float64, vz::Float64,
    speed::Float64, charge::Int, particle_rate_s1::Float64, ds::Float64,
)
    particle_rate_s1 > 0 || return
    mx, my, mz = 0.5(x + xnew), 0.5(y + ynew), 0.5(z + znew)
    cell = spatial_grid_cell(grid, mx, my, mz)
    isnothing(cell) && return
    ilon, ilat, ialt, volume = cell
    radius = sqrt(mx*mx + my*my + mz*mz)
    vr = (vx*mx + vy*my + vz*mz) / radius
    dt = ds / speed
    rate_over_volume = particle_rate_s1 / volume
    ic = charge + 1
    lock(grid.longitude_locks[ilon]) do
        grid.number_density_m3[ilon, ilat, ialt, ic] +=
            rate_over_volume * dt
        grid.total_flux_m2_s[ilon, ilat, ialt, ic] +=
            rate_over_volume * ds
        grid.signed_radial_flux_m2_s[ilon, ilat, ialt, ic] +=
            rate_over_volume * vr * dt
        grid.upward_radial_flux_m2_s[ilon, ilat, ialt, ic] +=
            rate_over_volume * max(vr, 0.0) * dt
        grid.downward_radial_flux_m2_s[ilon, ilat, ialt, ic] +=
            rate_over_volume * max(-vr, 0.0) * dt
    end
end

"""Accumulate one realized collision as a physical volume reaction rate."""
@inline function accumulate_spatial_reaction!(
    grid::SpatialGridAccumulator,
    x::Float64, y::Float64, z::Float64,
    charge_before::Int, target::Int, reaction::Int,
    energy_loss_ev::Float64, particle_rate_s1::Float64,
)
    particle_rate_s1 > 0 || return
    cell = spatial_grid_cell(grid, x, y, z)
    isnothing(cell) && return
    ilon, ilat, ialt, volume = cell
    rate = particle_rate_s1 / volume
    ic = charge_before + 1
    lock(grid.longitude_locks[ilon]) do
        grid.reaction_rate_m3_s1[
            ilon, ilat, ialt, ic, target, reaction
        ] += rate
        grid.reaction_event_count[ilon, ilat, ialt, reaction] += UInt32(1)
        grid.collision_energy_transfer_w_m3[ilon, ilat, ialt, ic] +=
            rate * energy_loss_ev * QE
    end
end

"""Deposit the remaining sub-cutoff projectile energy in its final cell."""
@inline function accumulate_spatial_thermalization!(
    grid::SpatialGridAccumulator,
    x::Float64, y::Float64, z::Float64,
    charge::Int, energy_ev::Float64, particle_rate_s1::Float64,
)
    particle_rate_s1 > 0 || return
    cell = spatial_grid_cell(grid, x, y, z)
    isnothing(cell) && return
    ilon, ilat, ialt, volume = cell
    lock(grid.longitude_locks[ilon]) do
        grid.cutoff_thermalization_w_m3[
            ilon, ilat, ialt, charge + 1
        ] += particle_rate_s1 * energy_ev * QE / volume
    end
end

"""Write grid geometry, moments, reactions, and energy to separate MAT files."""
function write_spatial_grid_mats(
    output_prefix::AbstractString,
    grid::SpatialGridAccumulator;
    metadata::AbstractDict=Dict{String,Any}(),
)
    prefix = abspath(output_prefix)
    mkpath(dirname(prefix))
    common = Dict{String,Any}(
        "format_version" => "marsaspen_spatial_grid_v1",
        "dimension_order" => "longitude, latitude, altitude, optional components",
        "longitude_edges_deg" => grid.longitude_edges_deg,
        "latitude_edges_deg" => grid.latitude_edges_deg,
        "altitude_edges_km" => grid.altitude_edges_km,
        "longitude_centers_deg" =>
            0.5 .* (grid.longitude_edges_deg[1:end-1] .+
                    grid.longitude_edges_deg[2:end]),
        "latitude_centers_deg" =>
            0.5 .* (grid.latitude_edges_deg[1:end-1] .+
                    grid.latitude_edges_deg[2:end]),
        "altitude_centers_km" =>
            0.5 .* (grid.altitude_edges_km[1:end-1] .+
                    grid.altitude_edges_km[2:end]),
        "coordinate_system" => "MSO",
    )
    merge!(common, Dict{String,Any}(string(k) => v for (k, v) in metadata))
    geometry_file = prefix * "_grid.mat"
    moments_file = prefix * "_moments.mat"
    reactions_file = prefix * "_reactions.mat"
    energy_file = prefix * "_energy.mat"
    matwrite(geometry_file, merge(copy(common), Dict(
        "cell_volume_m3" => grid.cell_volume_m3,
    )); compress=true)
    density_total = dropdims(sum(grid.number_density_m3; dims=4); dims=4)
    total_flux = dropdims(sum(grid.total_flux_m2_s; dims=4); dims=4)
    radial_flux = dropdims(
        sum(grid.signed_radial_flux_m2_s; dims=4); dims=4,
    )
    upward_flux = dropdims(
        sum(grid.upward_radial_flux_m2_s; dims=4); dims=4,
    )
    downward_flux = dropdims(
        sum(grid.downward_radial_flux_m2_s; dims=4); dims=4,
    )
    matwrite(moments_file, merge(copy(common), Dict(
        "charge_state_names" => ["H_ENA", "Hplus"],
        "number_density_by_charge_m3" => Float32.(grid.number_density_m3),
        "total_number_density_m3" => Float32.(density_total),
        "total_flux_by_charge_m2_s" => Float32.(grid.total_flux_m2_s),
        "total_flux_m2_s" => Float32.(total_flux),
        "signed_radial_flux_by_charge_m2_s" =>
            Float32.(grid.signed_radial_flux_m2_s),
        "signed_radial_flux_m2_s" => Float32.(radial_flux),
        "upward_radial_flux_by_charge_m2_s" =>
            Float32.(grid.upward_radial_flux_m2_s),
        "upward_radial_flux_m2_s" => Float32.(upward_flux),
        "downward_radial_flux_by_charge_m2_s" =>
            Float32.(grid.downward_radial_flux_m2_s),
        "downward_radial_flux_m2_s" => Float32.(downward_flux),
    )); compress=true)
    reaction_by_channel = dropdims(
        sum(grid.reaction_rate_m3_s1; dims=(4, 5)); dims=(4, 5),
    )
    ionization_by_target = dropdims(
        sum(grid.reaction_rate_m3_s1[:, :, :, :, :, 2]; dims=4);
        dims=4,
    )
    total_lya = dropdims(
        sum(grid.reaction_rate_m3_s1[:, :, :, :, :, 3]; dims=(4, 5));
        dims=(4, 5),
    )
    matwrite(reactions_file, merge(copy(common), Dict(
        "charge_state_names" => ["H_ENA", "Hplus"],
        "target_names" => ["CO2", "O", "N2"],
        "reaction_names" =>
            ["state_change", "ionization", "Ly_alpha", "elastic"],
        "reaction_rate_m3_s1" => Float32.(grid.reaction_rate_m3_s1),
        "reaction_rate_by_channel_m3_s1" =>
            Float32.(reaction_by_channel),
        "ionization_rate_by_target_m3_s1" =>
            Float32.(ionization_by_target),
        "total_lya_volume_emission_rate_photons_m3_s1" =>
            Float32.(total_lya),
        "raw_monte_carlo_event_count" => grid.reaction_event_count,
    )); compress=true)
    collision_energy = dropdims(
        sum(grid.collision_energy_transfer_w_m3; dims=4); dims=4,
    )
    cutoff_energy = dropdims(
        sum(grid.cutoff_thermalization_w_m3; dims=4); dims=4,
    )
    matwrite(energy_file, merge(copy(common), Dict(
        "charge_state_names" => ["H_ENA", "Hplus"],
        "collision_energy_transfer_by_charge_w_m3" =>
            Float32.(grid.collision_energy_transfer_w_m3),
        "collision_energy_transfer_w_m3" => Float32.(collision_energy),
        "cutoff_thermalization_by_charge_w_m3" =>
            Float32.(grid.cutoff_thermalization_w_m3),
        "cutoff_thermalization_w_m3" => Float32.(cutoff_energy),
        "total_energy_transfer_w_m3" =>
            Float32.(collision_energy .+ cutoff_energy),
    )); compress=true)
    (
        grid=geometry_file, moments=moments_file,
        reactions=reactions_file, energy=energy_file,
    )
end

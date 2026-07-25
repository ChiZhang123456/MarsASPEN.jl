module MarsASPEN

# Package entry point. Physics and numerical algorithms live in the included
# files below. The include order is significant because later files use types
# and helper functions defined by earlier files.
using MAT
using Random
using Statistics
using Base.Threads

export AspenModel, MonteCarloConfig, MonteCarloWeight, ParticleSummary,
       HistoryEvent, load_model,
       neutral_density, neutral_density_xyz, cartesian_to_lon_lat_alt,
       available_atmosphere_cases, run_ensemble, run_detailed_ensemble,
       run_binned_ensemble, run_phase_space_ensemble,
       run_density_crossing_ensemble,
       run_directional_flux_ensemble, run_radial_flux_ensemble,
       run_target_ionization_rate_ensemble,
       run_lya_volume_emission_ensemble,
       run_surface_diagnostic_map_ensemble,
       sample_injection_ensemble,
       run_spatial_grid_ensemble, write_spatial_grid_mats,
       thermal_speed_from_temperature_ev,
       maxwellian_importance_weight_3d, particle_density_weight,
       write_detailed_mat

include("types.jl")
include("monte_carlo_weight.jl")
include("cross_sections.jl")
include("initialization.jl")
include("atmosphere.jl")
include("spatial_grid.jl")
include("transport.jl")
include("ensembles.jl")
include("io.jl")

end

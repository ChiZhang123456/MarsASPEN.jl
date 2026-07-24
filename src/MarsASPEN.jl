module MarsASPEN

using MAT
using Random
using Statistics
using Base.Threads

export AspenModel, MonteCarloConfig, ParticleSummary, HistoryEvent, load_model,
       neutral_density, neutral_density_xyz, cartesian_to_lon_lat_alt,
       available_atmosphere_cases, run_ensemble, run_detailed_ensemble,
       run_binned_ensemble, run_phase_space_ensemble, write_detailed_mat

include("types.jl")
include("cross_sections.jl")
include("initialization.jl")
include("atmosphere.jl")
include("transport.jl")
include("ensembles.jl")
include("io.jl")

end

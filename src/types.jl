const QE = 1.602176634e-19
const AMU = 1.66053906660e-27
const MARS_RADIUS_KM = 3388.25
const H_MASS = 1.00782503223 * AMU
const HP_MASS = 1.007276466621 * AMU
const TARGET_MASS = (44.0095 * AMU, 15.999 * AMU, 28.0134 * AMU)
const ATMOSPHERE_MASS = (
    44.01 * AMU, 15.999 * AMU, 31.998 * AMU, 28.014 * AMU, 28.010 * AMU,
)
const ATMOSPHERE_SPECIES = (:CO2, :O, :O2, :N2, :CO)
const KB = 1.380649e-23
const MARS_G0 = 3.71
const REACTION_NAMES = (:state_change, :ionization, :lya, :elastic)
const MODEL_MIN_ALTITUDE_KM = 80.0

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
    min_altitude_km::Float64 = MODEL_MIN_ALTITUDE_KM
    max_altitude_km::Float64 = 1000.0
    max_collisions::Int = 2000
    max_steps_per_collision::Int = 100_000
    include_hot_o::Bool = true
    particle_weight::Float64 = 1.0
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

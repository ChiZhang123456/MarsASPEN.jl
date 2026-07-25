# Maxwellian importance sampling and physical macro-particle weights.
#
# The physical source distribution is f(v; U, T). Particles may be sampled
# from a broader distribution fs(v; U, Ts) to improve statistics in the tails.
# Each sampled particle then carries dimensionless importance weight f / fs.
# Physical density weights are Wn_i = n_source * W_i / sum(W), so their sum is
# exactly the prescribed source number density.

"""
Settings for Monte Carlo importance sampling and physical particle weights.

This object is intentionally separate from `MonteCarloConfig`. The config
controls particle transport, while this object controls how sampled
trajectories represent the physical source population.
"""
Base.@kwdef struct MonteCarloWeight
    # Sampling temperature divided by the physical source temperature.
    sampling_temperature_factor::Float64 = 1.0
    # Physical source number density. Zero selects unitless macro-particles.
    source_number_density_m3::Float64 = 0.0
    # Weight used when no physical source density is supplied.
    unit_particle_weight::Float64 = 1.0
end

"""Return Maxwellian thermal speed `sqrt(2 kT / m)` in m/s."""
@inline function thermal_speed_from_temperature_ev(
    temperature_ev::Real, mass_kg::Real=H_MASS,
)
    temperature_ev > 0 || throw(ArgumentError("temperature_ev must be positive"))
    mass_kg > 0 || throw(ArgumentError("mass_kg must be positive"))
    sqrt(2 * Float64(temperature_ev) * QE / Float64(mass_kg))
end

"""
Return the 3D Maxwellian importance weight `f / fs`.

`temperature_ev` is the physical source temperature and
`sampled_temperature_ev` is the temperature used to generate the velocity.
Density cancels in this ratio.
"""
@inline function maxwellian_importance_weight_3d(
    bulk_velocity_m_s::NTuple{3,<:Real},
    particle_velocity_m_s::NTuple{3,<:Real},
    temperature_ev::Real,
    sampled_temperature_ev::Real;
    mass_kg::Real=H_MASS,
)
    vth = thermal_speed_from_temperature_ev(temperature_ev, mass_kg)
    vsample = thermal_speed_from_temperature_ev(sampled_temperature_ev, mass_kg)
    dv2 = sum(
        (Float64(particle_velocity_m_s[i]) - Float64(bulk_velocity_m_s[i]))^2
        for i in 1:3
    )
    exp(3 * log(vsample / vth) + dv2 / vsample^2 - dv2 / vth^2)
end

"""
Convert one dimensionless importance weight to a density weight in m^-3.

`total_importance_weight` must be the sum over all source particles in the
ensemble. Summing the returned density weights therefore recovers
`source_density_m3`.
"""
@inline function particle_density_weight(
    importance_weight::Real,
    source_density_m3::Real,
    total_importance_weight::Real,
)
    importance_weight >= 0 ||
        throw(ArgumentError("importance_weight must be non-negative"))
    source_density_m3 >= 0 ||
        throw(ArgumentError("source_density_m3 must be non-negative"))
    total_importance_weight > 0 ||
        throw(ArgumentError("total_importance_weight must be positive"))
    Float64(source_density_m3) * Float64(importance_weight) /
    Float64(total_importance_weight)
end

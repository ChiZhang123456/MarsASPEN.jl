using MarsASPEN
using Test
using Statistics

const REPO = normpath(joinpath(@__DIR__, ".."))
const ATMOSPHERE_DIR = get(
    ENV,
    "MARSASPEN_ATMOSPHERE_DIR",
    joinpath(REPO, "data", "atmosphere"),
)
const MODEL = load_model(REPO; atmosphere_data_dir=ATMOSPHERE_DIR)

@testset "MGITM lower-atmosphere extrapolation" begin
    z1, z2 = MODEL.atmosphere.alt[1:2]
    for lon in (0.0, 45.0, 181.0), lat in (-45.0, 0.0, 45.0)
        n1 = MarsASPEN.density3(MODEL.atmosphere, lon, lat, z1)
        n2 = MarsASPEN.density3(MODEL.atmosphere, lon, lat, z2)
        n80 = MarsASPEN.density3(MODEL.atmosphere, lon, lat, 80.0)
        weight = (80.0 - z1) / (z2 - z1)
        expected = ntuple(
            i -> exp(log(n1[i]) + weight * (log(n2[i]) - log(n1[i]))), 5,
        )
        @test all(isapprox.(n80, expected; rtol=2e-13))
        @test MarsASPEN.density3(MODEL.atmosphere, lon, lat, 70.0) == n80

        t1 = MarsASPEN.temperature3(MODEL.atmosphere, lon, lat, z1)
        t2 = MarsASPEN.temperature3(MODEL.atmosphere, lon, lat, z2)
        expected_t80 = t1 + weight * (t2 - t1)
        @test MarsASPEN.temperature3(MODEL.atmosphere, lon, lat, 80.0) ≈
              expected_t80 rtol=2e-13
    end
    @test MonteCarloConfig().min_altitude_km == 80.0
end

@testset "Uniform dayside injection geometry" begin
    cfg = MonteCarloConfig(
        n_particles=10_000,
        seed=61,
        injection_geometry=:dayside_uniform,
        initial_charge_state=1,
        initial_temperature_ev=10.0,
    )
    weighting = MonteCarloWeight(source_number_density_m3=5.0e6)
    sample = sample_injection_ensemble(cfg; weighting=weighting)
    radius_km = sqrt.(sum(sample.position_m .^ 2; dims=2))[:] ./ 1000
    @test all(sample.position_m[:, 1] .>= 0)
    @test all(isapprox.(
        radius_km, MarsASPEN.MARS_RADIUS_KM + 600.0; atol=1e-9,
    ))
    @test minimum(sample.solar_zenith_angle_deg) >= 0
    @test maximum(sample.solar_zenith_angle_deg) <= 90
    @test abs(mean(sample.position_m[:, 1] ./ (radius_km .* 1000)) - 0.5) <
          0.01
    @test sum(sample.density_weight_m3) ≈ 5.0e6
    @test abs(mean(sample.velocity_m_s[:, 1]) / 1000 + 400) < 1.0
    @test abs(mean(sample.velocity_m_s[:, 2]) / 1000) < 1.0
    @test abs(mean(sample.velocity_m_s[:, 3]) / 1000) < 1.0
    @test !hasproperty(sample, :total_speed_flux_weight_m2_s)
end

@testset "Three-dimensional spatial diagnostics" begin
    cfg = MonteCarloConfig(
        n_particles=8,
        seed=62,
        injection_geometry=:dayside_uniform,
        initial_charge_state=1,
        initial_temperature_ev=10.0,
        min_altitude_km=590.0,
        max_step_m=1000.0,
    )
    result = run_spatial_grid_ensemble(
        MODEL, cfg;
        weighting=MonteCarloWeight(source_number_density_m3=5.0e6),
        altitude_edges_km=collect(590.0:1.0:600.0),
    )
    grid = result.grid
    @test size(grid.number_density_m3) == (72, 36, 10, 2)
    @test size(grid.reaction_rate_m3_s1) == (72, 36, 10, 2, 3, 4)
    @test all(grid.number_density_m3 .>= 0)
    @test all(grid.total_flux_m2_s .>= 0)
    @test all(grid.upward_radial_flux_m2_s .>= 0)
    @test all(grid.downward_radial_flux_m2_s .>= 0)
    @test grid.signed_radial_flux_m2_s ≈
          grid.upward_radial_flux_m2_s .-
          grid.downward_radial_flux_m2_s
    @test all(
        grid.total_flux_m2_s .+ 1e-12 .>=
        grid.upward_radial_flux_m2_s .+
        grid.downward_radial_flux_m2_s
    )
    @test sum(grid.number_density_m3) > 0
    @test sum(result.stop_counts) == cfg.n_particles
end

@testset "Local radial flux from density weights" begin
    cfg = MonteCarloConfig(
        n_particles=128,
        initial_charge_state=1,
        initial_temperature_ev=0.0,
        initial_speed_m_s=400_000.0,
        seed=31,
        include_hot_o=false,
    )
    weighting = MonteCarloWeight(source_number_density_m3=5.0e6)
    surfaces = collect(100.5:50.0:599.5)
    radial = run_radial_flux_ensemble(
        MODEL, cfg;
        weighting=weighting,
        altitude_surfaces_km=surfaces,
    )
    @test size(radial.radial_flux_m2_s) == (length(surfaces), 2, 2)
    @test all(radial.radial_flux_m2_s .>= 0)
    @test radial.signed_outward_flux_m2_s ≈
          radial.radial_flux_m2_s[:, :, 2] .-
          radial.radial_flux_m2_s[:, :, 1]
    @test radial.net_downward_flux_m2_s ≈
          -radial.signed_outward_flux_m2_s
    # At 599.5 km, the monoenergetic source is still H+ and travels radially
    # inward, so sum(Wn*|Vr|) equals n_source * 400 km/s.
    @test radial.radial_flux_m2_s[end, 2, 1] ≈ 2.0e12 rtol=1e-12

    doubled = run_radial_flux_ensemble(
        MODEL, cfg;
        weighting=MonteCarloWeight(source_number_density_m3=1.0e7),
        altitude_surfaces_km=surfaces,
    )
    @test doubled.radial_flux_m2_s ≈ 2 .* radial.radial_flux_m2_s
end

@testset "O ionization crossing estimator" begin
    cfg = MonteCarloConfig(
        n_particles=64,
        initial_charge_state=1,
        initial_temperature_ev=0.0,
        initial_speed_m_s=400_000.0,
        seed=33,
        include_hot_o=true,
    )
    surfaces = [200.5, 400.5, 599.5]
    result = run_target_ionization_rate_ensemble(
        MODEL, cfg;
        weighting=MonteCarloWeight(source_number_density_m3=5.0e6),
        altitude_surfaces_km=surfaces,
        target=:O,
    )
    @test size(result.ionization_rate_m3_s1) == (3, 2, 2)
    @test all(result.ionization_rate_m3_s1 .>= 0)
    @test result.ionization_rate_by_charge_m3_s1 ≈
          dropdims(sum(result.ionization_rate_m3_s1; dims=3); dims=3)
    @test result.total_ionization_rate_m3_s1 ≈
          vec(sum(result.ionization_rate_m3_s1; dims=(2, 3)))
    @test result.target_density_m3 ≈ [
        neutral_density(MODEL, 0.0, 0.0, altitude; include_hot_o=true).O
        for altitude in surfaces
    ]

    doubled = run_target_ionization_rate_ensemble(
        MODEL, cfg;
        weighting=MonteCarloWeight(source_number_density_m3=1.0e7),
        altitude_surfaces_km=surfaces,
        target=:O,
    )
    @test doubled.ionization_rate_m3_s1 ≈
          2 .* result.ionization_rate_m3_s1
end

@testset "Ly-alpha volume emission estimator" begin
    cfg = MonteCarloConfig(
        n_particles=32,
        initial_charge_state=0,
        initial_temperature_ev=0.0,
        initial_speed_m_s=400_000.0,
        seed=34,
        include_hot_o=true,
    )
    surfaces = [200.5, 400.5, 599.5]
    result = run_lya_volume_emission_ensemble(
        MODEL, cfg;
        weighting=MonteCarloWeight(source_number_density_m3=5.0e6),
        altitude_surfaces_km=surfaces,
    )
    @test size(result.volume_emission_rate_photons_m3_s1) == (3, 2, 3, 2)
    @test all(result.volume_emission_rate_photons_m3_s1 .>= 0)
    @test result.total_volume_emission_rate_photons_m3_s1 ≈
          vec(sum(result.volume_emission_rate_photons_m3_s1; dims=(2, 3, 4)))
    @test result.total_radiative_energy_rate_w_m3 ≈
          result.total_volume_emission_rate_photons_m3_s1 .*
          result.photon_energy_j
    @test result.photon_energy_j ≈
          6.62607015e-34 * 299_792_458.0 / 121.567e-9 rtol=1e-15
end

@testset "Atmosphere interpolation primitives" begin
    for altitude in (120.0, 200.0, 600.0)
        actual = MarsASPEN.density3(MODEL.atmosphere, 0.0, 0.0, altitude)
        @test all(isfinite, actual)
        @test all(>(0), actual)
    end
    for charge in (0, 1), target in 1:3
        @test MarsASPEN.sigma_at(MODEL.cross_sections, charge, target, 4, 1000.0) ≈
              7.848859363597747e-20 rtol=2e-14
    end
    for altitude in (100.0, 200.0, 600.0, 1000.0)
        value = MarsASPEN.hot_o_density(
            MODEL.hot_atmosphere, 0.0, 0.0, altitude,
        )
        @test isfinite(value)
        @test value > 0
    end
end

@testset "Complete neutral atmosphere interface" begin
    for solar in (70, 130, 200)
        model = load_model(REPO; solar=solar, ls=0)
        rho = neutral_density(model, 0.0, 0.0, 200.0)
        @test all(isfinite, (rho.CO2, rho.O, rho.O2, rho.N2, rho.CO, rho.Tn))
        @test rho.CO2 > 0
        @test rho.O > 0
        @test rho.Tn > 0
        # The processed GITM source files contain CO2 and O only.
        @test rho.O2 < 1e-250
        @test rho.N2 < 1e-250
        @test rho.CO < 1e-250
        rho_xyz = neutral_density_xyz(
            model, (3388.25 + 200.0) * 1000, 0.0, 0.0; position_unit=:m,
        )
        @test rho_xyz.CO2 ≈ rho.CO2 rtol=2e-13
    end
    low = neutral_density(load_model(REPO; solar=70, ls=0), 0.0, 0.0, 200.0)
    moderate = neutral_density(load_model(REPO; solar=130, ls=0), 0.0, 0.0, 200.0)
    high = neutral_density(load_model(REPO; solar=200, ls=0), 0.0, 0.0, 200.0)
    weight = (130 - 70) / (200 - 70)
    @test log(moderate.CO2) ≈
          (1 - weight) * log(low.CO2) + weight * log(high.CO2) rtol=2e-13
    @test moderate.Tn ≈
          (1 - weight) * low.Tn + weight * high.Tn rtol=2e-13
    @test length(available_atmosphere_cases()) == 12
    for case in available_atmosphere_cases()
        model = load_model(REPO; solar=case.f107, ls=case.ls)
        rho = neutral_density(model, 27.5, 0.0, 200.0)
        @test all(isfinite, (rho.CO2, rho.O, rho.O2, rho.N2, rho.CO, rho.Tn))
        @test all(>(0), (rho.CO2, rho.O, rho.O2, rho.N2, rho.CO, rho.Tn))
    end
end

@testset "Determinism and ensemble invariants" begin
    cfg = MonteCarloConfig(n_particles=20, seed=21)
    a = run_ensemble(MODEL, cfg; threaded=false)
    b = run_ensemble(MODEL, cfg; threaded=false)
    @test a == b
    @test all(r -> isfinite(r.final_energy_ev) && r.final_energy_ev >= 0, a)
    @test all(r -> r.n_collisions ==
        r.n_elastic + r.n_ionization + r.n_lya + r.n_state_change, a)
    binned = run_binned_ensemble(MODEL, cfg; altitude_edges_km=collect(80.0:10.0:1000.0))
    @test sum(binned.reaction_counts) == sum(r.n_collisions for r in binned.summaries)
end

@testset "Altitude-energy path-length histogram" begin
    cfg = MonteCarloConfig(n_particles=20, seed=17)
    altitude_edges = collect(80.0:10.0:600.0)
    energy_edges = 10.0 .^ range(1.0, 3.0, length=21)
    phase = run_phase_space_ensemble(
        MODEL, cfg;
        altitude_edges_km=altitude_edges,
        energy_edges_ev=energy_edges,
    )
    @test size(phase.path_length_m) ==
          (length(altitude_edges)-1, length(energy_edges)-1, 2)
    @test all(phase.path_length_m .>= 0)
    @test sum(phase.path_length_m[:,:,1]) > 0
    @test sum(phase.path_length_m[:,:,2]) > 0
    weighted = run_phase_space_ensemble(
        MODEL, MonteCarloConfig(n_particles=20, seed=17);
        weighting=MonteCarloWeight(unit_particle_weight=3.0),
        altitude_edges_km=altitude_edges,
        energy_edges_ev=energy_edges,
    )
    @test weighted.path_length_m ≈ 3 .* phase.path_length_m
end

@testset "Solar-wind proton directional flux" begin
    altitude_surfaces = collect(100.0:10.0:300.0)
    energy_edges = collect(10.0:20.0:2010.0)
    cfg = MonteCarloConfig(
        n_particles=40,
        seed=31,
        initial_charge_state=1,
        initial_temperature_ev=10.0,
    )
    weighting = MonteCarloWeight(source_number_density_m3=5.0e6)
    flux = run_directional_flux_ensemble(
        MODEL, cfg;
        weighting=weighting,
        altitude_surfaces_km=altitude_surfaces,
        energy_edges_ev=energy_edges,
    )
    @test size(flux.flux_m2_s) ==
          (length(altitude_surfaces), length(energy_edges)-1, 2, 2)
    @test all(flux.flux_m2_s .>= 0)
    @test sum(flux.flux_m2_s[:,:,2,1]) > 0
    @test sum(flux.flux_m2_s[end,:,2,1]) >
          sum(flux.flux_m2_s[end,:,1,1])
    @test sum(flux.stop_counts) == cfg.n_particles

    doubled = run_directional_flux_ensemble(
        MODEL, MonteCarloConfig(
            n_particles=40,
            seed=31,
            initial_charge_state=1,
            initial_temperature_ev=10.0,
        );
        weighting=MonteCarloWeight(source_number_density_m3=1.0e7),
        altitude_surfaces_km=altitude_surfaces,
        energy_edges_ev=energy_edges,
    )
    @test doubled.flux_m2_s ≈ 2 .* flux.flux_m2_s
end

@testset "Maxwellian importance and physical particle weights" begin
    @test !hasfield(MonteCarloConfig, :sampling_temperature_factor)
    @test !hasfield(MonteCarloConfig, :source_number_density_m3)
    @test !hasfield(MonteCarloConfig, :unit_particle_weight)
    default_weighting = MonteCarloWeight()
    @test default_weighting.sampling_temperature_factor == 1.0
    @test default_weighting.source_number_density_m3 == 0.0
    @test default_weighting.unit_particle_weight == 1.0

    bulk = (-400_000.0, 0.0, 0.0)
    velocity = (-380_000.0, 12_000.0, -8_000.0)
    @test maxwellian_importance_weight_3d(
        bulk, velocity, 10.0, 10.0; mass_kg=MarsASPEN.HP_MASS,
    ) ≈ 1.0

    vth = thermal_speed_from_temperature_ev(10.0, MarsASPEN.HP_MASS)
    vsample = thermal_speed_from_temperature_ev(50.0, MarsASPEN.HP_MASS)
    dv2 = sum((velocity[i] - bulk[i])^2 for i in 1:3)
    expected = (vsample / vth)^3 *
               exp(dv2 / vsample^2 - dv2 / vth^2)
    actual = maxwellian_importance_weight_3d(
        bulk, velocity, 10.0, 50.0; mass_kg=MarsASPEN.HP_MASS,
    )
    @test actual ≈ expected rtol=2e-14

    weights = [0.2, 0.7, 1.1, 2.0]
    density_weights = particle_density_weight.(
        weights, 5.0e6, sum(weights),
    )
    @test sum(density_weights) ≈ 5.0e6

    cfg = MonteCarloConfig(
        n_particles=40,
        seed=31,
        initial_charge_state=1,
        initial_temperature_ev=10.0,
    )
    weighting = MonteCarloWeight(
        sampling_temperature_factor=5.0,
        source_number_density_m3=5.0e6,
    )
    flux = run_directional_flux_ensemble(
        MODEL, cfg;
        weighting=weighting,
        altitude_surfaces_km=collect(100.0:20.0:300.0),
        energy_edges_ev=collect(10.0:40.0:2010.0),
    )
    @test flux.total_importance_weight > 0
    @test all(isfinite, flux.flux_m2_s)
    @test sum(flux.flux_m2_s) > 0

    density = run_density_crossing_ensemble(
        MODEL, cfg;
        weighting=weighting,
        altitude_surfaces_km=collect(100.5:20.0:280.5),
        energy_edges_ev=10.0 .^ range(0.0, 4.0, length=101),
    )
    @test size(density.density_weight_sum_m3) == (10, 100, 2)
    @test all(isfinite, density.density_weight_sum_m3)
    @test all(density.density_weight_sum_m3 .>= 0)
    @test sum(density.density_weight_sum_m3) > 0
end

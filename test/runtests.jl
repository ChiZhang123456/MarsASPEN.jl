using MarsASPEN
using Test

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

@testset "Python primitive parity" begin
    expected_density = (
        120.0 => (1.2113975782367725e17, 1.9307346587029375e15, 4.212254955721069e15),
        200.0 => (2.9749173651721215e13, 3.727716605690465e13, 1.380512797061012e13),
        600.0 => (0.018203028530616193, 1.0715626310704878e8, 2806.044609851673),
    )
    for (altitude, expected) in expected_density
        actual = MarsASPEN.density3(MODEL.atmosphere, 0.0, 0.0, altitude)
        @test all(isapprox.((actual[1], actual[2], actual[4]), expected; rtol=2e-13))
    end
    for charge in (0, 1), target in 1:3
        @test MarsASPEN.sigma_at(MODEL.cross_sections, charge, target, 4, 1000.0) ≈
              7.848859363597747e-20 rtol=2e-14
    end
    expected_hot_o = (
        100.0 => 9.806616016742045e8,
        200.0 => 2.144164213697179e9,
        600.0 => 7.446408390996631e8,
        1000.0 => 3.871697963971681e8,
    )
    for (altitude, expected) in expected_hot_o
        @test MarsASPEN.hot_o_density(MODEL.hot_atmosphere, 0.0, 0.0, altitude) ≈
              expected rtol=2e-13
    end
end

@testset "Complete neutral atmosphere interface" begin
    expected = Dict(
        70 => (
            2.9749173651721215e13, 3.727931022111834e13, 2.9396504567373694e11,
            1.380512797061012e13, 6.88442273662187e12, 194.91500000000002,
        ),
        130 => (
            1.1856414049182955e14, 6.914716189193995e13, 6.777959515736403e11,
            2.7057073784856965e13, 2.362370379565253e13, 267.1625,
        ),
        200 => (
            2.2218918625745278e14, 8.940117230789566e13, 9.23033917676682e11,
            3.3612170324763715e13, 4.46477869910226e13, 329.70125,
        ),
    )
    for solar in (70, 130, 200)
        model = load_model(REPO; solar=solar, ls=0)
        rho = neutral_density(model, 0.0, 0.0, 200.0)
        actual = (rho.CO2, rho.O, rho.O2, rho.N2, rho.CO, rho.Tn)
        @test all(isapprox.(actual, expected[solar]; rtol=2e-13))
        rho_xyz = neutral_density_xyz(
            model, (3388.25 + 200.0) * 1000, 0.0, 0.0; position_unit=:m,
        )
        @test rho_xyz.CO2 ≈ rho.CO2 rtol=2e-13
    end
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

using MarsASPEN
using Test

const REPO = normpath(joinpath(@__DIR__, ".."))
const ATMOSPHERE_DIR = get(
    ENV,
    "MARSASPEN_ATMOSPHERE_DIR",
    joinpath(REPO, "..", "py_aspen_github_package", "py_aspen", "neutral_density_model", "data"),
)
const MODEL = load_model(REPO; atmosphere_data_dir=ATMOSPHERE_DIR)

@testset "Python primitive parity" begin
    expected_density = (
        120.0 => (1.2113975782367725e17, 1.9307346587029375e15, 4.212254955721069e15),
        200.0 => (2.9749173651721215e13, 3.727716605690465e13, 1.380512797061012e13),
        600.0 => (0.018203028530616193, 1.0715626310704878e8, 2806.044609851673),
    )
    for (altitude, expected) in expected_density
        actual = MarsASPEN.density3(MODEL.atmosphere, 0.0, 0.0, altitude)
        @test all(isapprox.(actual, expected; rtol=2e-13))
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

@testset "Determinism and ensemble invariants" begin
    cfg = MonteCarloConfig(n_particles=20, seed=21)
    a = run_ensemble(MODEL, cfg; threaded=false)
    b = run_ensemble(MODEL, cfg; threaded=false)
    @test a == b
    @test all(r -> isfinite(r.final_energy_ev) && r.final_energy_ev >= 0, a)
    @test all(r -> r.n_collisions ==
        r.n_elastic + r.n_ionization + r.n_lya + r.n_state_change, a)
end

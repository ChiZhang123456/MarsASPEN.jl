# MarsASPEN examples

This folder contains two single-particle diagnostics and one weighted
100,000-particle Monte Carlo example. Run all commands from the repository root.

## 1. Single neutral H ENA

This case follows one monoenergetic 400 km/s H ENA from 600 km. It saves the
complete trajectory, energy, charge state, velocity, and collision history.
No physical source density or `MonteCarloWeight` object is needed.

```powershell
julia --project=. examples/run_h_ena_trajectory.jl
C:\Users\Win\.conda\envs\mars\python.exe examples/plot_single_trajectory.py examples/output/single_h_ena_400kms.mat
```

## 2. Single H+

This case uses the same altitude and speed but starts as H+. It is useful for
seeing charge exchange from H+ to H ENA and subsequent state changes.
No physical source density or `MonteCarloWeight` object is needed.

```powershell
julia --project=. examples/run_hplus_trajectory.jl
C:\Users\Win\.conda\envs\mars\python.exe examples/plot_single_trajectory.py examples/output/single_hplus_400kms.mat
```

Both single-particle figures contain altitude with reaction markers, energy,
charge state, and speed versus time. Reaction colors are shared between the
altitude and energy panels.

## 3. Weighted Monte Carlo example with 100,000 H ENA particles

This example injects 100,000 H ENA particles at 600 km with bulk velocity
`[-400, 0, 0]` km/s, physical temperature 10 eV, and source density 1 cm^-3.
It samples a broader 50 eV distribution and applies the importance weight
`f/fs` to every particle. The altitude-energy output uses 30 logarithmic
energy bins from 1 to 10,000 eV. The lower bound is 1 eV because a logarithmic
axis cannot include zero.

```powershell
julia --project=. -t auto examples/run_h_ena_100000_monte_carlo.jl
C:\Users\Win\.conda\envs\mars\python.exe examples/plot_h_ena_100000_monte_carlo.py examples/output/h_ena_100000_monte_carlo.mat
```

The resulting figure shows density-weighted track length versus altitude and
energy, separately for neutral H ENA and H+ created by electron stripping.

To inspect the H+ energy distribution at 550 km separately:

```powershell
C:\Users\Win\.conda\envs\mars\python.exe examples/plot_hplus_energy_at_altitude.py examples/output/h_ena_100000_monte_carlo.mat --altitude 550
```

The diagnostic shows both the native 5 eV bins and a less noisy 25 eV
rebinned curve.

The Julia code keeps transport settings and source weights separate:

```julia
config = MonteCarloConfig(
    n_particles=100_000,
    initial_speed_m_s=400_000.0,
    initial_charge_state=0,
    initial_temperature_ev=10.0,
)
weighting = MonteCarloWeight(
    sampling_temperature_factor=5.0,
    source_number_density_m3=1.0e6,
)

result = run_phase_space_ensemble(
    model, config; weighting=weighting,
)
```

When the sampling temperature differs from the physical temperature, each
velocity receives dimensionless importance weight

```text
W_i = f(v_i; U, T) / f_s(v_i; U, T_sample).
```

The physical density weight follows the py_aspen convention:

```text
Wn_i = n_source W_i / sum(W).
```

For flux through the injection plane, MarsASPEN uses

```text
Wflux_i = Wn_i max(-v_x,i, 0).
```

`MonteCarloWeight(sampling_temperature_factor=1)` samples the physical
Maxwellian directly, so all importance weights are one. Values above one
oversample velocity tails.

## Why `test/` is retained

The `test/` directory is not an example directory. It automatically verifies
the atmosphere interpolation, cross sections, deterministic random streams,
charge-state accounting, flux scaling, and Monte Carlo weight normalization.
It should remain in the package so future refactoring cannot silently change
the physics.

# MarsASPEN.jl

MarsASPEN.jl is a multithreaded Monte Carlo transport model for precipitating
H and H+ in the Martian atmosphere. It follows ASPEN-style collision physics
with CO2, O, and N2 and records energy deposition, ionization, charge-state
changes, H Ly-alpha production, and elastic scattering.

The repository contains:

* a typed Julia transport kernel under `src/`
* all 12 MGITM cold-atmosphere cases and all 12 MAMPS hot-O cases
* H and H+ collision cross-section tables under `data/cross_sections/`
* compact high-particle-count runs and detailed MAT v7.3 history output
* the `marsaspen-analysis` Python package under `analysis/`

## Source layout

The Julia package is separated by responsibility:

* `src/MarsASPEN.jl`: module entry point, exports, and include order
* `src/types.jl`: physical constants, model types, particle configuration, and outputs
* `src/initialization.jl`: MGITM, MAMPS, and cross-section model initialization
* `src/atmosphere.jl`: neutral-atmosphere interpolation and extrapolation
* `src/cross_sections.jl`: cross-section interpolation and reaction selection
* `src/monte_carlo_weight.jl`: Maxwellian importance sampling and physical particle weights
* `src/transport.jl`: one-particle propagation, collision kinematics, and state changes
* `src/ensembles.jl`: threaded ensembles and low-memory diagnostic histograms
* `src/spatial_grid.jl`: three-dimensional MSO lon-lat-alt diagnostics and multi-file MAT output
* `src/io.jl`: detailed MAT output

Simulation controls and source weights are deliberately separate:

```julia
config = MonteCarloConfig(
    n_particles=100_000,
    injection_geometry=:dayside_uniform,
    initial_charge_state=1,
    initial_speed_m_s=400_000.0,
    initial_temperature_ev=10.0,
)
weighting = MonteCarloWeight(
    sampling_temperature_factor=5.0,
    source_number_density_m3=5.0e6,
)
result = run_spatial_grid_ensemble(
    model, config; weighting=weighting,
)
```

`MonteCarloConfig` controls trajectories and numerical stopping conditions.
`MonteCarloWeight` controls importance sampling, source density normalization,
and the weight represented by each macro-particle.

A detailed Chinese guide to every source file, data convention, run script, and
analysis program is available in [`docs/CODE_GUIDE_ZH.md`](docs/CODE_GUIDE_ZH.md).

## Single-particle examples

The [`examples/`](examples/) directory contains reproducible 400 km/s H ENA and
H+ trajectories. These single-particle cases save detailed state, energy,
velocity, and collision history, then the shared Python program marks reaction
locations:

```powershell
julia --project=. examples/run_h_ena_trajectory.jl
julia --project=. examples/run_hplus_trajectory.jl
C:\Users\Win\.conda\envs\mars\python.exe examples/plot_single_trajectory.py examples/output/single_h_ena_400kms.mat
C:\Users\Win\.conda\envs\mars\python.exe examples/plot_single_trajectory.py examples/output/single_hplus_400kms.mat
```

Large-particle examples use a uniform 600 km MSO dayside hemispherical source,
not a single injection longitude and latitude:

```powershell
julia --project=. -t auto examples/sample_dayside_injection_100000.jl
julia --project=. -t auto examples/run_dayside_3d_100000.jl
C:\Users\Win\.conda\envs\mars\python.exe examples/plot_dayside_3d_maps_120km.py
C:\Users\Win\.conda\envs\mars\python.exe examples/plot_dayside_3d_altitude_profiles.py
C:\Users\Win\.conda\envs\mars\python.exe examples/plot_dayside_3d_sza_altitude.py
```

The complete source normalization, three-dimensional estimators, MAT contents,
and plotting commands are documented in
[`examples/README.md`](examples/README.md).

## Neutral atmosphere

The package includes MGITM and MAMPS data for Ls = 0, 90, 180, and 270 degrees
at F10.7 = 70, 130, and 200. MGITM supplies CO2, O, O2, N2, CO, and neutral
temperature. MAMPS supplies hot O. Above the MGITM top altitude, densities are
hydrostatically extrapolated using the local top-layer neutral temperature.
The native MGITM grid begins at 98.75 km. From 98.75 down to the 80 km transport
boundary, every longitude-latitude column is extrapolated linearly in log
density using its lowest two MGITM layers. Neutral temperature is extrapolated
linearly over the same interval. MAMPS hot O begins at 100 km and is zero below
its native lower boundary.

```julia
using MarsASPEN

model = load_model(solar="solar_min", ls=0)
rho = neutral_density(model, 0.0, 0.0, 200.0)
rho_xyz = neutral_density_xyz(
    model, (3388.25 + 200.0) * 1000, 0.0, 0.0; position_unit=:m
)
```

Returned densities are in m^-3 and `Tn` is in K. `O` is the sum of `O_cold`
from MGITM and `O_hot` from MAMPS.

## Julia setup and tests

```powershell
julia --project=. -e "using Pkg; Pkg.instantiate(); Pkg.test()"
```

## Uniform-dayside three-dimensional output

The production example injects 100,000 H+ macro particles uniformly in
spherical area over the 600 km MSO dayside hemisphere. Velocities are sampled
from a drifting Maxwellian with density $5~\mathrm{cm}^{-3}$, temperature
$10~\mathrm{eV}$, and MSO bulk velocity
$(-400,0,0)~\mathrm{km\,s}^{-1}$.

The simulation accumulates longitude, latitude, and 1 km altitude cells and
writes separate grid, moment, reaction, and energy MAT v7.3 files:

```powershell
julia --project=. -t auto examples/run_dayside_3d_100000.jl
```

The output supports longitude-latitude maps, altitude profiles, and
SZA-altitude maps of H and H+ density, scalar flux, radial flux, target
ionization, H Ly-alpha production, and projectile energy transfer.

## Python analysis

```powershell
C:\Users\Win\.conda\envs\mars\python.exe -m pip install -e analysis
marsaspen-plot output/aspen_10p_detailed.mat
```

The analysis package provides ordinary MAT and MAT v7.3 readers. The
uniform-dayside plotting examples read the gridded moment, reaction, and energy
files directly.

All non-mathematical figure text uses Arial.

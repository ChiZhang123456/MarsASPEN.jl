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
* `src/transport.jl`: one-particle propagation, collision kinematics, and state changes
* `src/ensembles.jl`: threaded ensembles and low-memory diagnostic histograms
* `src/io.jl`: detailed MAT output

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

## Compact ensemble benchmark

```powershell
julia --project=. -t auto scripts/benchmark.jl 1000000
```

On the development workstation, 1,000,000 particles completed in 131.0 s
using 28 Julia threads, corresponding to about 7,631 particles/s.

## Detailed MAT output

```powershell
julia --project=. -t auto scripts/run_detailed.jl 10 output/aspen_10p_detailed.mat
```

Detailed output uses flat arrays plus particle offsets in one compressed MAT
v7.3 file. Compact output is recommended for large ensembles. Saving every
transport step for 1,000,000 particles can require hundreds of GB.

## Altitude-binned reaction counts

Large ensembles can accumulate reaction counts directly into altitude bins
without storing every collision event:

```powershell
julia --project=. -t auto scripts/run_reaction_altitude_counts.jl 1000000 output/reaction_altitude_counts_1000000p.csv 10
C:\Users\Win\.conda\envs\mars\python.exe analysis/scripts/plot_reaction_altitude_counts.py output/reaction_altitude_counts_1000000p.csv
```

The CSV separates charge-state change, ionization, Lyman-alpha production, and
elastic collisions. Counts are collision events, so one particle can contribute
many events.

## H ENA and H+ altitude-energy distributions

```powershell
julia --project=. -t auto scripts/run_phase_space_histogram.jl 1000000 output/phase_space_1000000p.mat 1 100 1
C:\Users\Win\.conda\envs\mars\python.exe analysis/scripts/plot_phase_space_histogram.py output/phase_space_1000000p.mat
```

The final argument is the initial macro-particle weight, with a default of one.
Each trajectory segment contributes `particle_weight * ds` to its altitude,
energy, and charge-state bin. This path-length weighting avoids bias from the
adaptive transport step. A physical incident flux can be applied by replacing
the unit weight with the corresponding macro-particle weight.

## Python analysis

```powershell
C:\Users\Win\.conda\envs\mars\python.exe -m pip install -e analysis
marsaspen-plot output/aspen_10p_detailed.mat
```

All non-mathematical figure text uses Arial.

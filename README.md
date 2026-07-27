# MarsASPEN.jl

MarsASPEN.jl is a multithreaded Julia Monte Carlo model for transport of
solar wind H+ and energetic neutral H in the Martian atmosphere. It follows
individual three dimensional trajectories in MSO coordinates, samples
collisions with CO2, O, and N2, updates projectile energy and charge state
after every collision, and accumulates physical density, flux, ionization,
H Ly-alpha emission, and energy deposition diagnostics.

The package is designed as a complete Julia replacement for the transport
part of `py_aspen`. A separate Python package reads the Julia MAT output and
produces publication quality diagnostic figures.

## Current physical configuration

The standard uniform dayside simulation uses:

| Quantity | Value |
|---|---:|
| Injection altitude | 600 km |
| Injection surface | Uniform in area over the MSO dayside hemisphere |
| Initial projectile | H+ |
| MSO bulk velocity | \((-400,0,0)\) km s\(^{-1}\) |
| Physical temperature | 10 eV |
| Physical source density | 5 cm\(^{-3}\) |
| Lower transport boundary | 80 km |
| Upper transport boundary | 600 km |
| Thermalization cutoff | 10 eV |
| Collision count limit | None |
| Standard spatial grid | \(5^\circ\times5^\circ\times1\) km |

`max_collisions=nothing` means that a trajectory is not stopped because of
its physical collision count. A separate integration step safeguard remains
active to catch numerical pathologies.

## Repository structure

```text
MarsASPEN.jl/
├── src/                         Julia transport package
├── data/atmosphere/             GITM cold atmosphere and AMPS hot O
├── data/cross_sections/         H and H+ collision data
├── examples/                    Reproducible simulations and plotting scripts
├── analysis/marsaspen_analysis/ Python MAT readers and trajectory analysis
├── docs/                        Detailed Chinese technical documentation
└── test/                        Julia regression and physics consistency tests
```

The Julia source is separated by responsibility:

| File | Responsibility |
|---|---|
| `src/MarsASPEN.jl` | Module entry point and public API |
| `src/types.jl` | Constants, model containers, configuration, output records |
| `src/initialization.jl` | GITM, AMPS, and cross section loading |
| `src/atmosphere.jl` | Neutral atmosphere interpolation and extrapolation |
| `src/cross_sections.jl` | Cross section interpolation and event selection |
| `src/monte_carlo_weight.jl` | Maxwellian importance sampling and \(W_n\) |
| `src/transport.jl` | Free flight, collision sampling, scattering, state update |
| `src/ensembles.jl` | Threaded ensembles and low memory diagnostics |
| `src/spatial_grid.jl` | Three dimensional estimators and MAT v7.3 output |
| `src/io.jl` | Detailed variable length trajectory output |

## Installation

MarsASPEN requires Julia 1.10 or later.

```powershell
git clone https://github.com/ChiZhang123456/MarsASPEN.jl.git
cd MarsASPEN.jl
julia --project=. -e "using Pkg; Pkg.instantiate()"
```

Run all Julia tests with:

```powershell
julia --project=. -t 4 -e "using Pkg; Pkg.test()"
```

Install the optional Python analysis package with:

```powershell
C:\Users\Win\.conda\envs\mars\python.exe -m pip install -e analysis
```

## Load a model

```julia
using MarsASPEN

model = load_model(; solar="solar_min", ls=0)
```

Supported atmosphere selections are:

```julia
available_atmosphere_cases()
```

This returns combinations of \(L_s=0^\circ,90^\circ,180^\circ,270^\circ\)
and F10.7 = 70, 130, or 200.

Evaluate the atmosphere using MSO longitude, latitude, and altitude:

```julia
neutral = neutral_density(model, 0.0, 0.0, 200.0)
```

The returned number densities are in m\(^{-3}\), and `Tn` is in K. The
current processed GITM files contain CO2, cold O, and neutral temperature.
The O2, N2, and CO fields are set to a negligible numerical floor because
they are absent from these source files. AMPS supplies hot O inside its native
altitude range.

## Single particle examples

The H and H+ examples save every transport step and collision event:

```powershell
julia --project=. examples/run_h_ena_trajectory.jl
julia --project=. examples/run_hplus_trajectory.jl
```

Plot the trajectories with:

```powershell
C:\Users\Win\.conda\envs\mars\python.exe examples/plot_single_trajectory.py examples/output/single_h_ena_400kms.mat --altitude-min 100 --altitude-max 200 --output examples/figures/single_h_ena_400kms.png

C:\Users\Win\.conda\envs\mars\python.exe examples/plot_single_trajectory.py examples/output/single_hplus_400kms.mat --altitude-min 100 --altitude-max 200 --output examples/figures/single_hplus_400kms.png
```

These figures show trajectory altitude, energy, charge state, reaction times,
reaction counts, and sampled scattering angles.

## Uniform dayside source

The source position is sampled uniformly in spherical area over the dayside
hemisphere. Its global drift velocity is always specified in MSO coordinates.
The local radial velocity is therefore evaluated separately for every macro
particle.

Inspect 10,000,000 sampled source particles without running transport:

```powershell
julia --project=. -t 16 examples/sample_dayside_injection_10000000.jl
C:\Users\Win\.conda\envs\mars\python.exe examples/plot_dayside_injection_10000000.py
```

The stored macro particle weight is a number density weight:

$$
W_{n,i}
=
n_{\mathrm{source}}
\frac{w_i}{\sum_{p=1}^{N_{\mathrm{MC}}}w_p},
\qquad
[W_{n,i}]=\mathrm{m}^{-3}.
$$

It does not contain a radial speed or total speed factor. Velocity enters only
when a flux or a stationary source crossing rate is evaluated.

## Ten million particle simulation

Run the complete three dimensional example with:

```powershell
julia --project=. -t 16 examples/run_dayside_3d_10000000.jl
```

The simulation writes four compressed MAT v7.3 files:

```text
examples/output/dayside_hplus_10000000_3d_grid.mat
examples/output/dayside_hplus_10000000_3d_moments.mat
examples/output/dayside_hplus_10000000_3d_reactions.mat
examples/output/dayside_hplus_10000000_3d_energy.mat
```

The files contain grid geometry, H and H+ density and flux, target resolved
reaction rates, H Ly-alpha volume emission rate, and energy deposition rate.
Large MAT results under `examples/output/` are not tracked by Git.

Generate the main figures with:

```powershell
C:\Users\Win\.conda\envs\mars\python.exe examples/plot_dayside_3d_maps_120km.py
C:\Users\Win\.conda\envs\mars\python.exe examples/plot_dayside_3d_altitude_profiles.py
C:\Users\Win\.conda\envs\mars\python.exe examples/plot_dayside_3d_sza_altitude.py
C:\Users\Win\.conda\envs\mars\python.exe examples/plot_dayside_density_100_600km.py
```

## Three dimensional estimators

For a macro particle launch rate \(\dot N_i\), grid cell volume
\(V_{\mathrm{cell}}\), residence time \(\Delta t\), path length
\(\Delta s\), and radial velocity \(V_r\), MarsASPEN accumulates:

$$
n
=
\frac{1}{V_{\mathrm{cell}}}
\sum_i \dot N_i\Delta t_i,
$$

$$
F_{\mathrm{total}}
=
\frac{1}{V_{\mathrm{cell}}}
\sum_i \dot N_i\Delta s_i,
$$

$$
F_r
=
\frac{1}{V_{\mathrm{cell}}}
\sum_i \dot N_i V_{r,i}\Delta t_i.
$$

Upward and downward radial fluxes use
\(\max(V_r,0)\) and \(\max(-V_r,0)\), respectively. A realized collision
contributes:

$$
q_{\mathrm{event}}
=
\frac{\dot N_i}{V_{\mathrm{cell}}}.
$$

Raw Monte Carlo event counts are also stored, but they are not physical volume
rates.

## Documentation

Detailed documentation is organized as follows:

* [Documentation index](docs/README.md)
* [Code structure and file guide](docs/CODE_GUIDE_ZH.md)
* [Physics, collision sampling, and estimators](docs/PHYSICS_AND_ESTIMATORS_ZH.md)
* [MAT output schema](docs/OUTPUT_SCHEMA_ZH.md)
* [Running, convergence, and validation](docs/RUNNING_AND_VALIDATION_ZH.md)
* [Complete examples and plotting commands](examples/README.md)
* [Atmosphere data](data/atmosphere/README.md)
* [Cross section data](data/cross_sections/README.md)

## Reproducibility and known limitations

Each particle uses an RNG derived from `(seed, particle_id)`, so results are
independent of thread scheduling. The same seed and configuration reproduce
the same particle random streams.

Current limitations include:

* The processed GITM files contain only CO2, O, and neutral temperature.
* N2 cross sections are packaged, but the current atmosphere sets N2 to a
  negligible floor.
* MAMPS hot O is used only within its native altitude range.
* Magnetic and electric fields are not included.
* Gravity does not directly accelerate test particles during a free flight.
* H Ly-alpha output is a local volume emission rate. A line of sight
  integration is required to obtain Rayleigh brightness.
* The steady three dimensional normalization assumes the sampled injection
  surface and source model described in the technical documentation.

All nonmathematical figure text uses Arial.

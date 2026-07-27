# MarsASPEN documentation

The detailed technical documentation is written in Chinese so that model
assumptions, units, estimators, output variables, and validation steps can be
checked directly against the Julia source.

## Document map

| Document | Main content |
|---|---|
| [CODE_GUIDE_ZH.md](CODE_GUIDE_ZH.md) | Responsibilities and call relationships of every Julia and Python source file |
| [PHYSICS_AND_ESTIMATORS_ZH.md](PHYSICS_AND_ESTIMATORS_ZH.md) | Source sampling, \(W_n\), collision probability, state changes, scattering, density, flux, ionization, Ly-alpha, and energy deposition |
| [OUTPUT_SCHEMA_ZH.md](OUTPUT_SCHEMA_ZH.md) | Detailed trajectory MAT and four-file three dimensional MAT schemas |
| [RUNNING_AND_VALIDATION_ZH.md](RUNNING_AND_VALIDATION_ZH.md) | Installation, small tests, production runs, stop counts, convergence, and QA |

Additional data and example documentation:

* [Examples and plotting commands](../examples/README.md)
* [Atmosphere data](../data/atmosphere/README.md)
* [Cross section data](../data/cross_sections/README.md)

## Recommended reading order

For a first run:

1. Read the repository [README](../README.md).
2. Follow [RUNNING_AND_VALIDATION_ZH.md](RUNNING_AND_VALIDATION_ZH.md).
3. Use the commands in [examples/README.md](../examples/README.md).

For physical review:

1. Read [PHYSICS_AND_ESTIMATORS_ZH.md](PHYSICS_AND_ESTIMATORS_ZH.md).
2. Check atmosphere assumptions in
   [data/atmosphere/README.md](../data/atmosphere/README.md).
3. Check reaction channels in
   [data/cross_sections/README.md](../data/cross_sections/README.md).

For output analysis:

1. Read [OUTPUT_SCHEMA_ZH.md](OUTPUT_SCHEMA_ZH.md).
2. Install the Python package under `analysis/`.
3. Use the plotting programs documented in
   [examples/README.md](../examples/README.md).

## Current standard production configuration

```text
source: uniform MSO dayside hemisphere
injection altitude: 600 km
projectile: H+
bulk velocity: (-400, 0, 0) km/s
temperature: 10 eV
source density: 5 cm^-3
lower boundary: 80 km
upper boundary: 600 km
energy cutoff: 10 eV
collision limit: none
grid: 5 deg longitude x 5 deg latitude x 1 km altitude
```

When any of these values changes, record the new configuration in the output
metadata and in the simulation log.

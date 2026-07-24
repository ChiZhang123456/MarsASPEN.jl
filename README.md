# MarsASPEN.jl

MarsASPEN.jl is a multithreaded Monte Carlo transport model for precipitating
H and H+ in the Martian atmosphere. It follows ASPEN-style collision physics
with CO2, O, and N2 and records energy deposition, ionization, charge-state
changes, H Ly-alpha production, and elastic scattering.

The repository contains:

* a typed Julia transport kernel under `src/`
* MGITM cold atmosphere and MAMPS hot O interpolation
* H and H+ collision cross-section tables under `data/cross_sections/`
* compact high-particle-count runs and detailed MAT v7.3 history output
* the `marsaspen-analysis` Python package under `analysis/`

## Atmosphere data

The MGITM and MAMPS MAT files are not committed because the complete set is
about 144 MB. Set the atmosphere directory before running:

```powershell
$env:MARSASPEN_ATMOSPHERE_DIR = "F:\path\to\marsaspen_atmosphere"
```

The directory should contain files such as:

```text
gitm_ls000_f070.mat
mamps_ls000_f070.mat
```

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

## Python analysis

```powershell
C:\Users\Win\.conda\envs\mars\python.exe -m pip install -e analysis
marsaspen-plot output/aspen_10p_detailed.mat
```

All non-mathematical figure text uses Arial.

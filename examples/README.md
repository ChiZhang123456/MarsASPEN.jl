# MarsASPEN examples

The examples now use a uniform dayside hemispherical source for every
large-particle Monte Carlo simulation. The former 100,000-particle examples
that injected all particles at one longitude and latitude have been removed.
Single-particle examples remain available because they are intended to explain
individual trajectories and collision physics, not to represent a global
solar-wind source.

All commands below are run from the MarsASPEN.jl repository root.

## 1. Single H ENA trajectory

This example follows one neutral H particle with an initial speed of
400 km s\(^{-1}\). It records altitude, charge state, energy, reaction type,
scattering angle, and cumulative collision probability.

```powershell
julia --project=. examples/run_h_ena_trajectory.jl
C:\Users\Win\.conda\envs\mars\python.exe examples/plot_single_trajectory.py examples/output/single_h_ena_400kms.mat
```

The PNG is stored at
`examples/figures/single_h_ena_400kms.png`.

## 2. Single H+ trajectory

This example follows one proton with an initial speed of 400 km s\(^{-1}\).
Charge exchange can convert H+ into H ENA, and electron stripping can convert
H ENA back into H+.

```powershell
julia --project=. examples/run_hplus_trajectory.jl
C:\Users\Win\.conda\envs\mars\python.exe examples/plot_single_trajectory.py examples/output/single_hplus_400kms.mat
```

The PNG is stored at
`examples/figures/single_hplus_400kms.png`.

## 3. Atmospheric density cases

The packaged atmosphere combines MGITM cold neutral densities with MAMPS hot
oxygen. MGITM is logarithmically extrapolated from its lowest native grid
levels to 80 km at every longitude and latitude. MAMPS is used only inside its
available altitude range.

```powershell
C:\Users\Win\.conda\envs\mars\python.exe examples/plot_atmosphere_cases.py
```

The figure compares the available solar longitude and F10.7 cases from 80 to
1000 km:

`examples/figures/gitm_mamps_density_cases.png`.

The corresponding longitude-latitude maps at 150 km are generated with:

```powershell
C:\Users\Win\.conda\envs\mars\python.exe examples/plot_atmosphere_maps_150km.py
```

The upper four rows show MGITM CO2 and the lower four rows show total O,
defined as MGITM cold O plus MAMPS hot O. Density is interpolated in logarithmic
space to 150 km on the native 5 degree longitude-latitude grid. Each species
uses one common logarithmic color scale across all seasons and solar-activity
conditions:

`examples/figures/gitm_mamps_co2_o_maps_150km.png`.

## 4. Collision cross sections

The collision tables distinguish neutral H and H+ projectiles and the CO2, O,
and N2 targets. The channels include elastic scattering, charge-state change,
target ionization, Ly-alpha emission, and Balmer-alpha emission where data are
available.

```powershell
C:\Users\Win\.conda\envs\mars\python.exe examples/plot_cross_sections.py
```

The PNG is stored at
`examples/figures/collision_cross_sections.png`.

Cross sections are converted from cm\(^{2}\) to m\(^{2}\) before they are
used in SI collision rates.

## 5. Collision probability for a 1000 eV projectile

The following example evaluates collision probability every 1 km for a
1000 eV projectile moving downward through each atmosphere case:

```powershell
C:\Users\Win\.conda\envs\mars\python.exe examples/plot_1000ev_collision_probability.py
```

The local collision coefficient is

```math
\alpha(z)
=
\sum_j n_j(z)\sum_k \sigma_{j,k}(E),
```

where \(n_j\) is in m\(^{-3}\), \(\sigma_{j,k}\) is in m\(^{2}\), and
\(\alpha\) is in m\(^{-1}\). For a path length \(\Delta s\),

```math
P_{\mathrm{local}}
=
1-\exp[-\alpha(z)\Delta s].
```

The cumulative optical depth and collision probability are

```math
\tau(z)=\int_z^{1000\,\mathrm{km}}\alpha(s)\,ds,
\qquad
P_{\mathrm{cum}}=1-\exp[-\tau(z)].
```

The two PNG files are:

* `examples/figures/collision_probability_1000ev_local.png`
* `examples/figures/collision_probability_1000ev_cumulative.png`

## 6. Uniform dayside injection at 600 km

### Source geometry

The new large-particle source covers the complete MSO dayside hemisphere at
600 km altitude. Set

```julia
injection_geometry=:dayside_uniform
```

in `MonteCarloConfig`. MarsASPEN samples

```math
\mu=\cos(\mathrm{SZA})\sim U(0,1),
\qquad
\varphi\sim U(0,2\pi),
```

and constructs the MSO position

```math
\boldsymbol{r}
=
r_{\mathrm{inj}}
\left(
\mu,\,
\sqrt{1-\mu^2}\cos\varphi,\,
\sqrt{1-\mu^2}\sin\varphi
\right).
```

Here

```math
r_{\mathrm{inj}}=R_{\mathrm{Mars}}+600~\mathrm{km}.
```

Uniform \(\mu\) and azimuth give a uniform probability per unit spherical
surface area. All sampled positions satisfy \(x\geq0\) and
\(0\leq\mathrm{SZA}\leq90^\circ\).

### Velocity distribution

The physical proton distribution is a drifting three-dimensional Maxwellian:

```math
\boldsymbol{U}
=
(-400,0,0)~\mathrm{km\,s^{-1}},
\qquad
T=10~\mathrm{eV},
\qquad
n_{\mathrm{sw}}=5~\mathrm{cm^{-3}}.
```

The bulk velocity is specified in the global MSO frame. It is not rotated to
the local radial direction. Therefore, each particle has a different local
radial velocity,

```math
V_{r,i}
=
\frac{\boldsymbol{r}_i\mathbin{\cdot}\boldsymbol{v}_i}
{|\boldsymbol{r}_i|}.
```

At the subsolar point, the bulk velocity is almost entirely inward and
radial. Near the terminator, the same MSO velocity is nearly tangent to the
600 km injection sphere.

### Macro-particle normalization

The uniform-dayside simulation represents a stationary particle injection
rate through the spherical source surface. The area of the dayside hemisphere
is

```math
A_{\mathrm{day}}=2\pi r_{\mathrm{inj}}^2.
```

For physical distribution \(f(\boldsymbol{v})\), sampling distribution
\(f_s(\boldsymbol{v})\), and importance ratio

```math
w_i
=
\frac{f(\boldsymbol{v}_i)}
{f_s(\boldsymbol{v}_i)},
```

the inward number flux represented by the sampled ensemble is normalized with
the local inward speed

```math
V_{\mathrm{in},i}
=
\max(-V_{r,i},0).
```

The physical particle rate represented by macro particle \(i\) is

```math
\dot{N}_i
=
A_{\mathrm{day}}\,
n_{\mathrm{sw}}\,
\frac{w_iV_{\mathrm{in},i}}
{\sum_{p=1}^{N_{\mathrm{MC}}}w_p}.
```

Units are

```text
A_day:       m^2
n_sw:        m^-3
V_in:        m s^-1
w_i:         dimensionless
Ndot_i:      s^-1
```

Thus each simulated trajectory is a macro particle carrying a physical
particle rate, not one real proton and not a fixed point-source density.

### Inspecting the injection distribution

The following commands sample 100,000 initial positions and velocities without
running transport:

```powershell
julia --project=. -t auto examples/sample_dayside_injection_100000.jl
C:\Users\Win\.conda\envs\mars\python.exe examples/plot_dayside_injection_100000.py examples/output/dayside_hplus_injection_100000.mat
```

The resulting figure contains a three-dimensional Mars and injection-position
view together with the MSO \(V_x\), \(V_y\), and \(V_z\) distributions:

`examples/figures/dayside_hplus_injection_100000_4panel.png`.

## 7. Uniform-dayside three-dimensional Monte Carlo simulation

Run the complete 100,000-particle H+ simulation with:

```powershell
julia --project=. -t auto examples/run_dayside_3d_100000.jl
```

The example uses:

```text
initial altitude:       600 km
source geometry:        uniform dayside spherical area
initial species:        H+
MSO bulk velocity:      (-400, 0, 0) km s^-1
physical temperature:   10 eV
physical density:       5 cm^-3
number of particles:    100,000
longitude bin width:    5 degrees
latitude bin width:     5 degrees
altitude bin width:     1 km
altitude range:         80 to 600 km
```

### Three-dimensional estimators

For grid-cell volume \(V_{\mathrm{cell}}\), trajectory residence time
\(\Delta t\), path length \(\Delta s=|\boldsymbol{v}|\Delta t\), and
macro-particle rate \(\dot N_i\), MarsASPEN accumulates

```math
n
=
\sum_i
\frac{\dot N_i\Delta t_i}{V_{\mathrm{cell}}},
```

```math
F_{\mathrm{total}}
=
\sum_i
\frac{\dot N_i\Delta s_i}{V_{\mathrm{cell}}},
```

```math
F_r
=
\sum_i
\frac{\dot N_iV_{r,i}\Delta t_i}{V_{\mathrm{cell}}}.
```

The corresponding units are m\(^{-3}\), m\(^{-2}\) s\(^{-1}\), and
m\(^{-2}\) s\(^{-1}\). Upward and downward radial fluxes are accumulated
separately using \(\max(V_r,0)\) and \(\max(-V_r,0)\).

Each realized reaction contributes

```math
q_{\mathrm{event}}
=
\frac{\dot N_i}{V_{\mathrm{cell}}}
```

in m\(^{-3}\) s\(^{-1}\). Reactions are stored separately by projectile
charge state, atmospheric target, and reaction channel.

### MAT output files

The simulation writes four MAT v7.3 files:

* `dayside_hplus_100000_3d_grid.mat` contains longitude, latitude, and
  altitude edges and centers, plus exact spherical cell volumes.
* `dayside_hplus_100000_3d_moments.mat` contains total and charge-resolved
  number density, total scalar flux, signed radial flux, upward radial flux,
  and downward radial flux.
* `dayside_hplus_100000_3d_reactions.mat` contains reaction rates and raw
  Monte Carlo event counts by charge, target, and channel. It also stores
  target-resolved ionization rates and total H Ly-alpha volume emission.
* `dayside_hplus_100000_3d_energy.mat` contains collision energy transfer,
  sub-10 eV thermalization, and their sum in W m\(^{-3}\).

Julia writes arrays in longitude, latitude, altitude order followed by any
component dimensions. Some HDF5 readers expose dimensions in reverse order.
The supplied Python scripts identify dimensions from coordinate lengths
instead of assuming NumPy axis order.

## 8. Uniform-dayside diagnostic figures

### Initial positions and velocities

```powershell
C:\Users\Win\.conda\envs\mars\python.exe examples/plot_dayside_injection_100000.py examples/output/dayside_hplus_injection_100000.mat
```

Output:

`examples/figures/dayside_hplus_injection_100000_4panel.png`

### Longitude and latitude maps at 120 km

```powershell
C:\Users\Win\.conda\envs\mars\python.exe examples/plot_dayside_3d_maps_120km.py
```

The eight panels show number density, total scalar flux, downward radial flux,
upward radial flux, signed radial flux, total target ionization rate,
H Ly-alpha volume emission rate, and projectile energy transfer rate:

`examples/figures/dayside_hplus_100000_3d_maps_120km_8panel.png`

### Altitude profiles

```powershell
C:\Users\Win\.conda\envs\mars\python.exe examples/plot_dayside_3d_altitude_profiles.py
```

Profiles cover 100 to 300 km and compare the global spherical-area mean with
the dayside-area mean:

`examples/figures/dayside_hplus_100000_3d_altitude_profiles_8panel.png`

Latitude weighting uses the exact spherical factor

```math
\sin(\phi_{\mathrm{upper}})
-
\sin(\phi_{\mathrm{lower}}).
```

### SZA and altitude distributions

```powershell
C:\Users\Win\.conda\envs\mars\python.exe examples/plot_dayside_3d_sza_altitude.py
```

The MSO solar zenith angle is

```math
\mathrm{SZA}
=
\cos^{-1}(\cos\phi\cos\lambda).
```

The script uses 5 degree SZA bins, 1 km altitude bins, and the 100 to 300 km
altitude range. Every SZA annulus is averaged with exact spherical cell-area
weights. The dotted line at SZA \(=90^\circ\) marks the terminator.

Output:

`examples/figures/dayside_hplus_100000_3d_sza_altitude_8panel.png`

## 9. Analysis package

The Python package under `analysis/marsaspen_analysis` provides MAT and MAT
v7.3 readers for trajectory and three-dimensional output. The former
point-source radial-flux, ionization-profile, and Ly-alpha-profile modules
have been removed. Uniform-dayside diagnostics are calculated directly from
the gridded moments, reactions, and energy MAT files.

Install the reader in editable mode with:

```powershell
C:\Users\Win\.conda\envs\mars\python.exe -m pip install -e analysis
```

## 10. Tests

The Julia tests verify atmosphere interpolation, cross sections, Monte Carlo
weighting, uniform-dayside injection geometry, and three-dimensional spatial
diagnostics:

```powershell
julia --project=. test/runtests.jl
```

The remaining Python tests verify ordinary MAT and MAT v7.3 input:

```powershell
C:\Users\Win\.conda\envs\mars\python.exe -m pytest analysis/tests
```

# MarsASPEN examples

This folder contains two single-particle diagnostics and two weighted
100,000-particle Monte Carlo examples. Run all commands from the repository
root. The source files contain detailed comments describing physical units,
charge-state conventions, importance sampling, histogram dimensions, MAT
fields, and limitations of the crossing-weight estimator.

## 1. Single neutral H ENA

This case follows one monoenergetic 400 km/s H ENA from 600 km. It saves the
complete trajectory, energy, charge state, velocity, and collision history.
No physical source density or `MonteCarloWeight` object is needed.

```powershell
julia --project=. examples/run_h_ena_trajectory.jl
C:\Users\Win\.conda\envs\mars\python.exe examples/plot_single_trajectory.py examples/output/single_h_ena_400kms.mat
```

The plotting script displays only the atmospheric segment from 80 to 400 km
and saves a 600 dpi PNG in `examples/figures`.

## 2. Single H+

This case uses the same altitude and speed but starts as H+. It is useful for
seeing charge exchange from H+ to H ENA and subsequent state changes.
No physical source density or `MonteCarloWeight` object is needed.

```powershell
julia --project=. examples/run_hplus_trajectory.jl
C:\Users\Win\.conda\envs\mars\python.exe examples/plot_single_trajectory.py examples/output/single_hplus_400kms.mat
```

The corresponding publication-style PNG is stored in `examples/figures`.

Both single-particle figures use a 2 by 2 quantitative layout containing
altitude with reaction markers, energy, charge state, and speed versus time.
Reaction colors are shared between altitude and energy. Each reference figure
is stored as a 600 dpi PNG.

## 3. Weighted Monte Carlo example with 100,000 H ENA particles

This example injects 100,000 numerical H ENA particles at 600 km with bulk
velocity `[-400, 0, 0]` km/s, physical temperature 10 eV, and source number
density

```math
n_{\mathrm{source}} = 1~\mathrm{cm^{-3}}
                    = 10^{6}~\mathrm{m^{-3}}.
```

Each simulated trajectory is a **macro particle**, not one individual
physical hydrogen atom. Macro particle $i$ carries a density weight
$W_{n,i}$, which specifies the share of the physical source number density
represented by that trajectory:

```math
[W_{n,i}] = \mathrm{m^{-3}}.
```

If a physical source volume $\Delta V$, in $\mathrm{m^3}$, is specified,
the number of real particles represented by macro particle $i$ in that
volume is

```math
N_{\mathrm{real},i} = W_{n,i}\Delta V,
```

which is dimensionless because it is a particle count. MarsASPEN stores
$W_{n,i}$, rather than assuming a particular source volume.

### Physical and sampling velocity distributions

The requested source is a drifting three-dimensional Maxwellian. Its
normalized velocity probability density is

```math
p(\boldsymbol{v};\boldsymbol{U},T)=
\left(\frac{m}{2\pi k_{\mathrm{B}}T}\right)^{3/2}
\exp\left[
-\frac{m|\boldsymbol{v}-\boldsymbol{U}|^2}
       {2k_{\mathrm{B}}T}
\right],
```

where

* $\boldsymbol{v}$ is particle velocity, in $\mathrm{m\,s^{-1}}$;
* $\boldsymbol{U}$ is bulk velocity, in $\mathrm{m\,s^{-1}}$;
* $T$ is expressed in the input as eV, with $k_{\mathrm{B}}T=T_{\mathrm{eV}}q_{\mathrm{e}}$, in J;
* $m$ is the hydrogen mass, in kg;
* $p(\boldsymbol{v})$ has units $(\mathrm{m\,s^{-1}})^{-3}=\mathrm{s^3\,m^{-3}}$;
* $\int p(\boldsymbol{v})\,d^3v=1$.

The physical phase-space number-density distribution is therefore

```math
f_n(\boldsymbol{v})=
n_{\mathrm{source}}p(\boldsymbol{v};\boldsymbol{U},T),
```

with units $\mathrm{m^{-3}}(\mathrm{m\,s^{-1}})^{-3}=\mathrm{s^3\,m^{-6}}$.

Sampling every numerical particle directly from the physical 10 eV
distribution is valid. However, this example deliberately samples from a
broader Maxwellian with

```math
T_{\mathrm{s}}=50~\mathrm{eV},
```

so that the velocity tails contain more numerical trajectories. Let the
normalized sampling probability density be

```math
p_{\mathrm{s}}(\boldsymbol{v})=
p(\boldsymbol{v};\boldsymbol{U},T_{\mathrm{s}}).
```

### Importance weight

For a sampled velocity $\boldsymbol{v}_i$, the dimensionless importance
weight is

```math
w_i=
\frac{p(\boldsymbol{v}_i;\boldsymbol{U},T)}
       {p_{\mathrm{s}}(\boldsymbol{v}_i;\boldsymbol{U},T_{\mathrm{s}})}.
```

For the two drifting Maxwellians used here, MarsASPEN evaluates this as

```math
w_i=
\left(\frac{v_{\mathrm{th,s}}}{v_{\mathrm{th}}}\right)^3
\exp\left[
\frac{|\boldsymbol{v}_i-\boldsymbol{U}|^2}{v_{\mathrm{th,s}}^2}
-
\frac{|\boldsymbol{v}_i-\boldsymbol{U}|^2}{v_{\mathrm{th}}^2}
\right],
```

where

```math
v_{\mathrm{th}}=\sqrt{\frac{2T_{\mathrm{eV}}q_{\mathrm{e}}}{m}},
\qquad
v_{\mathrm{th,s}}=
\sqrt{\frac{2T_{\mathrm{s,eV}}q_{\mathrm{e}}}{m}}.
```

Both probability densities have the same units, so

```math
[w_i]=1.
```

The importance weight corrects the deliberately broadened sampling
distribution back to the requested physical 10 eV distribution.

### Density weight carried by each macro particle

For $N_{\mathrm{MC}}$ simulated particles, MarsASPEN converts $w_i$ into
the physical density weight

```math
\boxed{
W_{n,i}=
n_{\mathrm{source}}
\frac{w_i}{\displaystyle\sum_{j=1}^{N_{\mathrm{MC}}}w_j}
}
```

with

```math
[W_{n,i}]=\mathrm{m^{-3}}.
```

This normalization guarantees

```math
\sum_{i=1}^{N_{\mathrm{MC}}}W_{n,i}=
n_{\mathrm{source}}.
```

Thus, the 100,000 numerical trajectories collectively represent the complete
physical source density. They do not represent only 100,000 real atoms. If
the sampling and physical temperatures are identical, all $w_i=1$, and the
formula reduces to

```math
W_{n,i}=\frac{n_{\mathrm{source}}}{N_{\mathrm{MC}}}.
```

For the present example, this equal-weight reference value would be

```math
\frac{10^6~\mathrm{m^{-3}}}{10^5}=
10~\mathrm{m^{-3}}
```

per macro particle. Because importance sampling is used, the actual
$W_{n,i}$ values are unequal, but their sum remains exactly
$10^6~\mathrm{m^{-3}}$.

### Using Wn in model diagnostics

For an altitude-energy histogram, every crossing of macro particle $i$ is
added to its corresponding altitude, energy, and charge-state bin using
$W_{n,i}$:

```math
H_{a,e,q}=
\sum_{i\in(a,e,q)}W_{n,i},
\qquad
[H_{a,e,q}]=\mathrm{m^{-3}}.
```

This example directly plots that accumulated density weight. It does not
multiply by trajectory path length and does not divide by energy-bin width.
Consequently, the plotted quantity is density weight accumulated per discrete
altitude-energy bin, not differential density per eV.

For a vertical number-flux diagnostic, the local radial velocity must also be
included:

```math
V_{r,i}=
\frac{\boldsymbol{r}_i\cdot\boldsymbol{v}_i}
       {|\boldsymbol{r}_i|},
\qquad
[V_{r,i}]=\mathrm{m\,s^{-1}},
```

```math
F_i=W_{n,i}|V_{r,i}|,
\qquad
[F_i]=\mathrm{m^{-2}\,s^{-1}}.
```

Downward and upward crossings are accumulated separately. The signed outward
flux is $F_{\mathrm{up}}-F_{\mathrm{down}}$.

### Running the example

The altitude-energy output uses 100 logarithmic energy bins from 1 to
10,000 eV. The lower bound is 1 eV because a logarithmic axis cannot include
zero. The displayed energy range is 10 to 3,000 eV.

```powershell
julia --project=. -t auto examples/run_h_ena_100000_monte_carlo.jl
C:\Users\Win\.conda\envs\mars\python.exe examples/plot_h_ena_100000_monte_carlo.py examples/output/h_ena_100000_monte_carlo.mat
```

The resulting figure separately shows neutral H ENA and H+ created by
electron stripping. A rendered reference figure is stored at
`examples/figures/h_ena_100000_altitude_energy.png`.

## 4. Weighted Monte Carlo example with 100,000 H+ particles

This case injects 100,000 solar-wind protons at 600 km with a bulk velocity of
`[-400, 0, 0]` km/s, physical temperature 10 eV, and physical density
5 cm^-3. The importance sampler again uses 50 eV, and the f/fs correction
restores the requested 10 eV distribution. Charge exchange produces the H ENA
population shown in the left panel.

```powershell
julia --project=. -t auto examples/run_hplus_100000_monte_carlo.jl
C:\Users\Win\.conda\envs\mars\python.exe examples/plot_h_ena_100000_monte_carlo.py examples/output/hplus_100000_monte_carlo.mat --output examples/figures/hplus_100000_altitude_energy.png
```

The plotting script reads `initial_species`, source density, temperature,
altitude, speed, and particle count from the MAT file, so the title is not
hard-coded for the neutral example. For the proton-source figure, the H ENA
color range is 1e1 to 1e6 m^-3 and the dominant H+ range is 1e1 to 1e7 m^-3.
A rendered reference figure is stored at
`examples/figures/hplus_100000_altitude_energy.png`.

## 5. MGITM and MAMPS density profiles

This Nature-style 4 by 3 panel figure compares all packaged atmosphere cases:
four seasons, $L_s=0,90,180,270$ degrees, and three F10.7 values, 70, 130,
and 200. Profiles are evaluated at 0 degrees longitude and 0 degrees latitude
from 80 to 1000 km. The five cold species use MGITM, including the same lower
log-density and upper hydrostatic extrapolations as the Julia model. Hot O uses
MAMPS only inside its native altitude range.

```powershell
C:\Users\Win\.conda\envs\mars\python.exe examples/plot_atmosphere_cases.py
```

The PNG is stored at
`examples/figures/gitm_mamps_density_cases.png`.

## 6. Collision cross sections

The cross-section figure has one row for H ENA and one row for H+, with columns
for CO2, O, and N2. Every panel shows state change, target ionization,
Ly-alpha production, and elastic scattering. Cross sections are converted from
the source-table unit cm^2 to the model SI unit m^2.

```powershell
C:\Users\Win\.conda\envs\mars\python.exe examples/plot_cross_sections.py
```

The PNG is stored at
`examples/figures/collision_cross_sections.png`.

## 7. Fixed 1000 eV collision probability every 1 km

This diagnostic holds the projectile energy fixed at 1000 eV and computes the
total collision coefficient for H ENA and H+ at every 1 km altitude:

```text
alpha(z) = sum_j n_j(z) sum_k sigma_jk(1000 eV)
P_local(z) = 1 - exp[-alpha(z) 1000 m]
P_cumulative(z) = 1 - exp[-integral_z^1000km alpha(s) ds]
```

The Python calculation reproduces Julia's longitude and latitude interpolation,
lower and upper atmosphere extrapolation, MAMPS range handling, and internal
768-point cross-section interpolation. It was checked directly against
`MarsASPEN.local_state` at 100, 200, 600, and 1000 km.

```powershell
C:\Users\Win\.conda\envs\mars\python.exe examples/plot_1000ev_collision_probability.py
```

The two PNG outputs are:

- `examples/figures/collision_probability_1000ev_local.png`
- `examples/figures/collision_probability_1000ev_cumulative.png`

Because energy is deliberately fixed, these curves describe first-collision
optical depth rather than a complete energy-degrading Monte Carlo trajectory.

## 8. Local radial flux profile

This example injects 100,000 H+ particles at 600 km with 400 km/s bulk speed,
10 eV temperature, and 5 cm^-3 density. At every one-kilometer spherical
altitude surface, it evaluates the local radial velocity and accumulates

```text
Vr_i = v_i dot r_hat
Wn_i = n_source (f/fs)_i / sum_j(f/fs)_j       [m^-3]
F_i = Wn_i abs(Vr_i)                           [m^-2 s^-1]
```

Downward and upward magnitudes are saved separately for H ENA and H+. The
signed outward radial flux is `F_upward - F_downward`, so negative values
indicate net precipitation and positive values indicate net escape.

```powershell
julia --project=. -t auto examples/run_hplus_100000_radial_flux.jl
C:\Users\Win\.conda\envs\mars\python.exe examples/plot_radial_flux_profile.py examples/output/hplus_100000_radial_flux.mat --output examples/figures/hplus_100000_radial_flux_profile.png
```

The PNG is stored at
`examples/figures/hplus_100000_radial_flux_profile.png`.

The equivalent initially neutral H ENA source uses the same density,
temperature, bulk speed, altitude grid, and flux definition:

```powershell
julia --project=. -t auto examples/run_h_ena_100000_radial_flux.jl
C:\Users\Win\.conda\envs\mars\python.exe examples/plot_radial_flux_profile.py examples/output/h_ena_100000_radial_flux.mat --output examples/figures/h_ena_100000_radial_flux_profile.png
```

The reusable Python calculations are in
`analysis/marsaspen_analysis/flux.py`. They provide local vertical velocity,
single-particle `Wn * Vr`, and directional and net altitude profiles.

To inspect the H+ energy distribution at 550 km separately:

```powershell
C:\Users\Win\.conda\envs\mars\python.exe examples/plot_hplus_energy_at_altitude.py examples/output/h_ena_100000_monte_carlo.mat --altitude 550
```

To compare the highest saved altitude layer with the analytical 10 eV
drifting Maxwellian injection distribution, run:

```powershell
C:\Users\Win\.conda\envs\mars\python.exe examples/plot_energy_at_600km.py examples/output/h_ena_100000_monte_carlo.mat
```

The highest row represents the 599.5 km crossing surface. Its current
density-weight histogram combines downward and upward crossings, so returned
particles can produce a low-energy tail even this close to the injection
boundary.

The one-altitude diagnostic shows all native logarithmic bins and can combine
neighboring bins with its `--rebin` option to reduce Monte Carlo noise.

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

result = run_density_crossing_ensemble(
    model, config;
    weighting=weighting,
    altitude_surfaces_km=altitude_centers_km,
    energy_edges_ev=energy_edges_ev,
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

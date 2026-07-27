# Collision cross sections and scattering data

This directory contains the collision data used by MarsASPEN for energetic
neutral hydrogen, H, and protons, H+, interacting with the main neutral
atmospheric targets CO2, O, and N2.

The transport code reads these files in
[`src/cross_sections.jl`](../../src/cross_sections.jl). The single-particle
collision update is implemented in
[`src/transport.jl`](../../src/transport.jl).
The complete probability and estimator definitions are documented in
[`docs/PHYSICS_AND_ESTIMATORS_ZH.md`](../../docs/PHYSICS_AND_ESTIMATORS_ZH.md).

## Files

| Projectile | Target | File |
|---|---|---|
| H | CO2 | `H_CO2_cross_sections.txt` |
| H | O | `H_O_cross_sections.txt` |
| H | N2 | `H_N2_cross_sections.txt` |
| H+ | CO2 | `Hplus_CO2_cross_sections.txt` |
| H+ | O | `Hplus_O_cross_sections.txt` |
| H+ | N2 | `Hplus_N2_cross_sections.txt` |
| H or H+ | atmospheric target | `scattering_angle_distribution.txt` |

Each cross-section table covers projectile energies from 1 eV to 9 MeV.
Energy is stored in eV and cross section is stored in cm$^2$ in the text files.
MarsASPEN converts cross sections to SI units during loading:

$$
1\ \mathrm{cm^2}=10^{-4}\ \mathrm{m^2}.
$$

## Neutral H tables

The H-impact files have the following columns:

```text
energy_eV sigma_el_cm2 sigma_01_cm2 sigma_La_cm2 sigma_ia_cm2 sigma_Ha_cm2
```

| Column | Meaning | MarsASPEN channel |
|---|---|---|
| `energy_eV` | H projectile kinetic energy | Energy coordinate |
| `sigma_el_cm2` | Elastic scattering | Elastic |
| `sigma_01_cm2` | Electron stripping, H becomes H+ | State change |
| `sigma_La_cm2` | H Ly-alpha production | Ly-alpha |
| `sigma_ia_cm2` | Ionization of the atmospheric target by H | Ionization |
| `sigma_Ha_cm2` | Balmer-alpha production | Stored in the table, not currently used |

For example, electron stripping changes the projectile charge state according
to

$$
\mathrm{H}+j\rightarrow\mathrm{H^+}+j+e^-,
$$

where $j$ is CO2, O, or N2.

## Proton tables

The H+-impact files have the following columns:

```text
energy_eV sigma_10_cm2 sigma_ip_cm2 sigma_La_cm2 sigma_el_cm2 sigma_Ha_cm2
```

| Column | Meaning | MarsASPEN channel |
|---|---|---|
| `energy_eV` | H+ projectile kinetic energy | Energy coordinate |
| `sigma_10_cm2` | Charge exchange, H+ becomes neutral H | State change |
| `sigma_ip_cm2` | Ionization of the atmospheric target by H+ | Ionization |
| `sigma_La_cm2` | H Ly-alpha production | Ly-alpha |
| `sigma_el_cm2` | Elastic scattering | Elastic |
| `sigma_Ha_cm2` | Balmer-alpha production | Stored in the table, not currently used |

Charge exchange is represented schematically as

$$
\mathrm{H^+}+j\rightarrow\mathrm{H}+j^+.
$$

## Common internal channel order

The H and H+ files use different column orders. MarsASPEN rearranges both into
one internal order:

| Internal reaction index | Reaction |
|---:|---|
| 1 | State change |
| 2 | Atmospheric target ionization |
| 3 | H Ly-alpha production |
| 4 | Elastic collision |

The internal cross-section array is

```text
sigma[projectile charge state, target, reaction, energy]
```

with:

- projectile charge state 1: neutral H
- projectile charge state 2: H+
- target 1: CO2
- target 2: O
- target 3: N2

Julia uses one-based array indices. Elsewhere in the particle records, the
physical charge-state value is 0 for H and 1 for H+.

## Energy interpolation

MarsASPEN constructs an internal logarithmic grid containing 768 points from
1 eV to 10 MeV. Each tabulated cross section is linearly interpolated in
energy onto this grid. No cross-section extrapolation is applied. A requested
energy outside the tabulated range returns zero.

The fixed inelastic energy losses are read from the `Q = ... eV` comments in
each file. MarsASPEN stores the absolute loss magnitude and subtracts it after
the reaction:

$$
E_{\mathrm{after}}
=
\max\left(E_{\mathrm{before}}-\lvert Q\rvert,0\right).
$$

Elastic energy transfer is calculated separately using two-body collision
kinematics and the sampled laboratory-frame scattering angle.

## Collision probability and event selection

At position $\mathbf r$ and projectile energy $E$, the total collision
coefficient is

$$
\alpha(\mathbf r,E)
=
\sum_j n_j(\mathbf r)
\sum_k \sigma_{j,k}(E),
$$

where:

- $n_j$ is target number density in m$^{-3}$
- $\sigma_{j,k}$ is the cross section in m$^2$
- $j$ identifies CO2, O, or N2
- $k$ identifies the four modeled reaction channels
- $\alpha$ has units m$^{-1}$

For a path segment of length $\Delta s$, the collision probability is

$$
P_{\mathrm{collision}}
=
1-\exp\left[-\alpha(\mathbf r,E)\Delta s\right].
$$

Once a collision occurs, the target and reaction channel are sampled using
their relative $n_j\sigma_{j,k}$ contributions:

$$
P(j,k\mid\mathrm{collision})
=
\frac{n_j\sigma_{j,k}}
{\displaystyle\sum_{j'}\sum_{k'}n_{j'}\sigma_{j',k'}}.
$$

The projectile energy and charge state are updated immediately. The collision
coefficient is then recalculated before the next free-flight segment.

## Scattering-angle distribution

`scattering_angle_distribution.txt` contains two columns:

```text
random_number scattering_angle_lab_deg
```

The first column is a cumulative random coordinate and the second column is
the corresponding laboratory-frame polar scattering angle. MarsASPEN draws

$$
u\sim U(0,1)
$$

and linearly interpolates the table to obtain

$$
\theta_{\mathrm{LAB}}=\theta(u).
$$

The azimuthal angle is sampled independently:

$$
\phi\sim U(0,2\pi).
$$

The pair $(\theta_{\mathrm{LAB}},\phi)$ rotates the velocity direction after
every realized collision. Therefore, scattering is stochastic even when the
projectile energy, target, and reaction channel are the same.

## Data provenance

The header of every text file records the source MATLAB figure from which the
curve was extracted. The cross-section curves originate from the Figure 3
panels, and the inverse-CDF-style scattering distribution originates from
Figure 2. These header comments are retained so each numerical table remains
traceable to its source curve.

## Current limitations

* The same laboratory-frame scattering-angle lookup table is currently used
  for all modeled targets, energies, reaction channels, and projectile charge
  states.
* Balmer-alpha columns are retained in the source tables but are not included
  in the four-channel transport kernel.
* N2 cross sections are available, but the packaged processed GITM atmosphere
  does not contain a physical N2 density field.
* Cross sections outside the tabulated energy range are zero.
* H Ly-alpha production is treated as a local event channel. Radiative
  transfer and resonant reabsorption are not included.

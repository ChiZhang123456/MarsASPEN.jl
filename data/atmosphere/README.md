# Neutral atmosphere data

This directory contains the processed GITM cold atmosphere and AMPS hot atomic
oxygen files read directly by MarsASPEN. No intermediate `_sph_grid.mat` layer
is used.

The transport model combines:

$$
n_{\mathrm O,\mathrm{total}}
=
n_{\mathrm O,\mathrm{GITM}}
+
n_{\mathrm O,\mathrm{AMPS}}.
$$

GITM provides CO2, cold O, and neutral temperature. AMPS provides the hot O
corona.

## Coordinates

The original files are already expressed in Mars Solar Orbital coordinates.
The subsolar point is

$$
\lambda_{\mathrm{MSO}}=0^\circ,\qquad
\phi_{\mathrm{MSO}}=0^\circ.
$$

No geographic-to-MSO rotation is applied by MarsASPEN. The GITM horizontal
grid uses cell centers at longitudes $2.5^\circ, 7.5^\circ, \ldots, 357.5^\circ$
and latitudes $-87.5^\circ, -82.5^\circ, \ldots, 87.5^\circ$. Consequently,
the closest GITM cell center to the mathematical subsolar point has a small,
nonzero SZA. The AMPS longitude grid includes $0^\circ$, but its latitude
cell centers are also offset by $2.5^\circ$.

For any MSO longitude $\lambda$ and latitude $\phi$, the solar zenith
angle used by the examples is:

$$
\mathrm{SZA}
=
\cos^{-1}\left(\cos\phi\cos\lambda\right).
$$

Longitude interpolation is periodic across $0^\circ$ and $360^\circ$.
Latitude interpolation is linear between neighboring cell centers.

## Case mapping

| Ls | GITM source | AMPS source |
|---:|---|---|
| $0^\circ$ | `gitm_aequ*_subL0_alt220.mat` | `dsmc_aequ*.mat` |
| $90^\circ$ | `gitm_aph*_subL0_alt220.mat` | `dsmc_aph*.mat` |
| $180^\circ$ | `gitm_aequ*_subL180_alt220.mat` | `dsmc_aequ*.mat` |
| $270^\circ$ | `gitm_per*_subL270_alt220.mat` | `dsmc_per*.mat` |

The suffix `min` supplies the F10.7 = 70 case and `max` supplies the
F10.7 = 200 case. Because the source collection has no independent moderate
solar-activity file, the F10.7 = 130 density fields are interpolated
logarithmically between min and max. Neutral temperature is interpolated
linearly.

The supported combinations returned by `available_atmosphere_cases()` are:

```text
Ls = 0, 90, 180, 270 degrees
F10.7 = 70, 130, 200
```

## File naming

GITM files follow:

```text
gitm_<season><activity>_<subsolar-tag>_alt220.mat
```

AMPS files follow:

```text
dsmc_<season><activity>.mat
```

Here:

| Token | Meaning |
|---|---|
| `aequ` | Equinox case |
| `aph` | Aphelion case |
| `per` | Perihelion case |
| `min` | F10.7 = 70 |
| `max` | F10.7 = 200 |
| `subL0`, `subL180`, `subL270` | Source subsolar longitude tag |

## Fields and units

GITM supplies `NCO2`, `NO`, and `temp`. Number densities are converted during
loading from $\mathrm{cm}^{-3}$ to $\mathrm{m}^{-3}$, and temperature is in K. The supplied
processed GITM files do not contain N2, O2, or CO, so these fields are set to
a negligible positive floor and do not contribute to transport collisions.

AMPS supplies hot atomic oxygen density as `dens_oh`. It is converted during
loading from $\mathrm{cm}^{-3}$ to $\mathrm{m}^{-3}$. MarsASPEN reads these original field
names directly and does not create another repackaged atmosphere data layer.

The coordinate variables read from both products are:

```text
Longitude
Latitude
Altitude
```

MarsASPEN reshapes the vector form into longitude, latitude, altitude cubes.
The duplicated periodic longitude column is removed, leaving 72 longitude
cells and 36 latitude cells.

## Solar activity interpolation

For a positive density $n$, F10.7 = 130 is evaluated in logarithmic space:

$$
\ln n_{130}
=
(1-w)\ln n_{70}+w\ln n_{200},
$$

where:

$$
w=\frac{130-70}{200-70}.
$$

Neutral temperature uses linear interpolation:

$$
T_{130}
=
(1-w)T_{70}+wT_{200}.
$$

## Spatial interpolation

Cold neutral densities are stored and interpolated in log density. At a point
inside the native grid, MarsASPEN performs periodic longitude interpolation,
linear latitude interpolation, and linear altitude interpolation in
$\ln n$.

Neutral temperature is interpolated linearly.

## Lower atmosphere extension

The native GITM grid begins near 98.75 km. For every longitude and latitude
column, MarsASPEN extrapolates from the lowest two native altitude layers down
to 80 km:

$$
\ln n(h)
=
\ln n(h_1)
+
\frac{h-h_1}{h_2-h_1}
\left[
\ln n(h_2)-\ln n(h_1)
\right].
$$

Temperature uses the analogous linear expression in $T$, rather than
$\ln T$. Requests below 80 km are clamped to the 80 km atmosphere because
80 km is the model lower boundary.

## Above the GITM top

Cold species above the GITM top use a hydrostatic exponential extension:

$$
n_j(h)
=
n_j(h_{\mathrm{top}})
\exp\left[
-\frac{h-h_{\mathrm{top}}}{H_j}
\right],
$$

with:

$$
H_j
=
\frac{k_BT_{\mathrm{top}}}{m_jg_{\mathrm{top}}}.
$$

The gravity used in the scale height is:

$$
g_{\mathrm{top}}
=
g_0
\left(
\frac{R_{\mathrm{Mars}}}
{R_{\mathrm{Mars}}+h_{\mathrm{top}}}
\right)^2.
$$

The standard transport upper boundary is 600 km, but the hydrostatic extension
also supports atmosphere diagnostic plots above that altitude.

## AMPS altitude behavior

AMPS hot O is interpolated only inside its native altitude range. It is set to
zero below the first AMPS altitude and above the last AMPS altitude. MarsASPEN
does not extrapolate the hot O corona beyond the supplied AMPS grid.

## Transport targets

The public `neutral_density` function returns:

```text
CO2, O, O2, N2, CO, Tn, O_cold, O_hot
```

The current processed GITM files do not contain O2, N2, or CO. These species
are assigned a negligible positive floor. The collision kernel currently
evaluates CO2, O, and N2 because cross section tables exist for those targets.
In the packaged atmosphere, the effective transport targets are therefore CO2
and total O.

## Loading examples

```julia
using MarsASPEN

model = load_model(; ls=0, solar="solar_min")
neutral = neutral_density(model, 0.0, 0.0, 150.0)

println(neutral.CO2)
println(neutral.O_cold)
println(neutral.O_hot)
println(neutral.O)
println(neutral.Tn)
```

An external atmosphere directory can be selected with:

```julia
model = load_model(
    ;
    ls=0,
    solar=70,
    atmosphere_data_dir="path/to/data/atmosphere",
)
```

Alternatively, set the environment variable:

```text
MARSASPEN_ATMOSPHERE_DIR
```

## Validation

After replacing any MAT file, verify:

1. Longitude and latitude arrays are monotonic.
2. Longitude covers one complete 360 degree period.
3. Density values are positive after unit conversion.
4. Temperature is positive.
5. The 80 km extension joins continuously to the native GITM grid.
6. The $0^\circ$ longitude and $0^\circ$ latitude point is the subsolar
   direction.
7. AMPS hot O is zero outside its native altitude range.

## Limitations

* The files are processed model products, not native GITM or AMPS restart
  outputs.
* F10.7 = 130 is interpolated, not independently simulated.
* O2, N2, and CO are unavailable in the packaged GITM products.
* No geographic to MSO time dependent rotation is applied.
* Atmospheric time variability during a particle trajectory is not included.

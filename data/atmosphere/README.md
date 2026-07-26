# Neutral atmosphere data

The `GITM/` and `AMPS/` subdirectories are direct copies of the MATLAB files
from:

```text
D:\Work_Work\Mars\MAVEN\test_particle_jl\neutral\GITM
D:\Work_Work\Mars\MAVEN\test_particle_jl\neutral\AMPS
```

## Coordinates

The original files are already expressed in Mars Solar Orbital coordinates.
The subsolar point is

$$
\lambda_{\mathrm{MSO}}=0^\circ,\qquad
\phi_{\mathrm{MSO}}=0^\circ.
$$

No geographic-to-MSO rotation is applied by MarsASPEN. The GITM horizontal
grid uses cell centers at longitudes \(2.5^\circ, 7.5^\circ,\ldots,357.5^\circ\)
and latitudes \(-87.5^\circ,-82.5^\circ,\ldots,87.5^\circ\). Consequently,
the closest GITM cell center to the mathematical subsolar point has a small,
nonzero SZA. The AMPS longitude grid includes \(0^\circ\), but its latitude
cell centers are also offset by \(2.5^\circ\).

## Case mapping

| Ls | GITM source | AMPS source |
|---:|---|---|
| \(0^\circ\) | `gitm_aequ*_subL0_alt220.mat` | `dsmc_aequ*.mat` |
| \(90^\circ\) | `gitm_aph*_subL0_alt220.mat` | `dsmc_aph*.mat` |
| \(180^\circ\) | `gitm_aequ*_subL180_alt220.mat` | `dsmc_aequ*.mat` |
| \(270^\circ\) | `gitm_per*_subL270_alt220.mat` | `dsmc_per*.mat` |

The suffix `min` supplies the F10.7 = 70 case and `max` supplies the
F10.7 = 200 case. Because the source collection has no independent moderate
solar-activity file, the F10.7 = 130 density fields are interpolated
logarithmically between min and max. Neutral temperature is interpolated
linearly.

## Fields and units

GITM supplies `NCO2`, `NO`, and `temp`. Number densities are converted during
loading from cm\(^{-3}\) to m\(^{-3}\), and temperature is in K. The supplied
processed GITM files do not contain N2, O2, or CO, so these fields are set to
a negligible positive floor and do not contribute to transport collisions.

AMPS supplies hot atomic oxygen density as `dens_oh`. It is converted during
loading from cm\(^{-3}\) to m\(^{-3}\). MarsASPEN reads these original field
names directly and does not create another repackaged atmosphere data layer.

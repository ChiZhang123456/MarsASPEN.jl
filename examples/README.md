# MarsASPEN examples

These examples run one 400 km/s projectile from 600 km and save the complete
trajectory, energy, charge state, and collision history.

## Neutral H ENA

```powershell
julia --project=. examples/run_h_ena_trajectory.jl
C:\Users\Win\.conda\envs\mars\python.exe examples/plot_single_trajectory.py examples/output/single_h_ena_400kms.mat
```

## H+

```powershell
julia --project=. examples/run_hplus_trajectory.jl
C:\Users\Win\.conda\envs\mars\python.exe examples/plot_single_trajectory.py examples/output/single_hplus_400kms.mat
```

The Python figure contains altitude with reaction markers, energy, charge
state, and speed versus time. Reaction colors are shared between the altitude
and energy panels.

## Monte Carlo source weights

Large ensembles use the routines in `src/monte_carlo_weight.jl`.

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

`sampling_temperature_factor = 1` samples the physical Maxwellian directly,
so all importance weights are one. Values above one oversample velocity tails.
The solar-wind production script uses a factor of five.

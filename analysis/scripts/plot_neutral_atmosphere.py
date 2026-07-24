"""Plot one MGITM cold-atmosphere column together with MAMPS hot O."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from scipy.constants import Boltzmann
from scipy.interpolate import RegularGridInterpolator
from scipy.io import loadmat

mpl.rcParams["font.family"] = "Arial"
mpl.rcParams["mathtext.fontset"] = "dejavusans"

MARS_RADIUS_KM = 3389.5
MARS_G0 = 3.72076
AMU_KG = 1.66053906660e-27
SPECIES = ("CO2", "O", "O2", "N2", "CO")
MASS_AMU = np.array((44.0, 16.0, 32.0, 28.0, 28.0))


def periodic_log_profile(
    values: np.ndarray,
    lon: np.ndarray,
    lat: np.ndarray,
    alt: np.ndarray,
    lon_query: float,
    lat_query: float,
    altitude: np.ndarray,
    log_values: bool = True,
    extrapolate_lower: bool = False,
) -> np.ndarray:
    lon_extended = np.r_[lon[-1] - 360.0, lon, lon[0] + 360.0]
    values_extended = np.concatenate(
        (values[-1:, :, :], values, values[:1, :, :]), axis=0
    )
    query_lon = ((lon_query - lon_extended[0]) % 360.0) + lon_extended[0]
    interpolator = RegularGridInterpolator(
        (lon_extended, lat, alt),
        np.log(np.maximum(values_extended, 1e-300)) if log_values else values_extended,
        bounds_error=False,
        fill_value=None,
    )
    points = np.column_stack(
        (
            np.full_like(altitude, query_lon),
            np.full_like(altitude, np.clip(lat_query, lat[0], lat[-1])),
            np.minimum(altitude, alt[-1])
            if extrapolate_lower
            else np.clip(altitude, alt[0], alt[-1]),
        )
    )
    result = interpolator(points)
    return np.exp(result) if log_values else result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--atmosphere-dir", type=Path, required=True)
    parser.add_argument("--ls", type=int, default=0)
    parser.add_argument("--f107", type=int, default=70)
    parser.add_argument("--lon", type=float, default=0.0)
    parser.add_argument("--lat", type=float, default=0.0)
    parser.add_argument("--altitude-range", type=float, nargs=2, default=(80, 1000))
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    tag = f"ls{args.ls:03d}_f{args.f107:03d}"
    gitm = loadmat(args.atmosphere_dir / f"gitm_{tag}.mat")
    mamps = loadmat(args.atmosphere_dir / f"mamps_{tag}.mat")
    z = np.arange(args.altitude_range[0], args.altitude_range[1] + 1.0, 1.0)

    lon_g = np.ravel(gitm["lon_deg"])
    lat_g = np.ravel(gitm["lat_deg"])
    alt_g = np.ravel(gitm["alt_km"])
    cold = {}
    for species in SPECIES:
        cold[species] = periodic_log_profile(
            np.asarray(gitm[f"n{species}"], dtype=float),
            lon_g, lat_g, alt_g, args.lon, args.lat, z,
            extrapolate_lower=True,
        )

    temperature = periodic_log_profile(
        np.asarray(gitm["Tn"], dtype=float),
        lon_g, lat_g, alt_g, args.lon, args.lat, z,
        log_values=False,
        extrapolate_lower=True,
    )
    above = z > alt_g[-1]
    if np.any(above):
        dz_km = z[above] - alt_g[-1]
        gravity = MARS_G0 * (MARS_RADIUS_KM / (MARS_RADIUS_KM + alt_g[-1])) ** 2
        for species, mass in zip(SPECIES, MASS_AMU):
            scale_height_km = (
                Boltzmann * temperature[above][0] / (mass * AMU_KG * gravity) / 1000
            )
            cold[species][above] *= np.exp(-dz_km / scale_height_km)
        temperature[above] = temperature[np.flatnonzero(~above)[-1]]

    hot_o = periodic_log_profile(
        np.asarray(mamps["nO_hot"], dtype=float),
        np.ravel(mamps["lon_deg"]),
        np.ravel(mamps["lat_deg"]),
        np.ravel(mamps["alt_km"]),
        args.lon,
        args.lat,
        z,
    )
    hot_alt = np.ravel(mamps["alt_km"])
    hot_o[(z < hot_alt[0]) | (z > hot_alt[-1])] = 0.0

    fig, axes = plt.subplots(1, 2, figsize=(11, 7), sharey=True, constrained_layout=True)
    colors = {
        "CO2": "#d62728", "O": "#1f77b4", "O2": "#9467bd",
        "N2": "#ff7f0e", "CO": "#8c564b",
    }
    for species in SPECIES:
        axes[0].plot(cold[species], z, lw=2, color=colors[species], label=species)
    axes[0].plot(
        np.where(hot_o > 0, hot_o, np.nan), z,
        lw=2.2, color="#2ca02c", label="Hot O, MAMPS",
    )
    axes[0].plot(
        cold["O"] + hot_o, z, lw=2, color="black", ls="--", label="Total O"
    )
    axes[0].set_xscale("log")
    axes[0].set_xlabel(r"Neutral density (m$^{-3}$)")
    axes[0].set_ylabel("Altitude (km)")
    axes[0].legend(loc="best")

    axes[1].plot(temperature, z, lw=2, color="#e377c2")
    axes[1].set_xlabel("MGITM neutral temperature (K)")
    for ax in axes:
        ax.axhline(
            alt_g[0], color="0.5", lw=1.3, ls="--",
            label="MGITM lower edge, extrapolation below" if ax is axes[1] else None,
        )
        ax.axhline(
            alt_g[-1], color="0.35", lw=1.3, ls=":",
            label="MGITM top, extrapolation above" if ax is axes[1] else None,
        )
        ax.set_ylim(*args.altitude_range)
        ax.grid(True, which="both", color="0.85", lw=0.7)
    axes[1].legend(loc="best")
    fig.suptitle(
        f"MGITM + MAMPS atmosphere, Ls = {args.ls} deg, "
        f"F10.7 = {args.f107}, lon = {args.lon:g} deg, lat = {args.lat:g} deg"
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=220)
    plt.close(fig)
    print(f"MGITM altitude grid: {alt_g[0]:g} to {alt_g[-1]:g} km")
    print(
        "MAMPS altitude grid: "
        f"{np.ravel(mamps['alt_km'])[0]:g} to {np.ravel(mamps['alt_km'])[-1]:g} km"
    )
    print(f"output={args.output.resolve()}")


if __name__ == "__main__":
    main()

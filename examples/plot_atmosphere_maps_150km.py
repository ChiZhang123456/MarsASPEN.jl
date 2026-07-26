"""Plot one native-longitude GITM and AMPS atmosphere case at 150 km."""

from __future__ import annotations

from pathlib import Path
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "analysis"))
from marsaspen_analysis import load_atmosphere_case  # noqa: E402

DATA = REPO / "data" / "atmosphere"
OUTPUT = (
    REPO / "examples" / "figures" /
    "gitm_mamps_co2_o_maps_150km.png"
)

LS_DEG = 0
F107 = 200
ALTITUDE_KM = 150.0

mpl.rcParams.update({
    "font.family": "Arial",
    "mathtext.fontset": "dejavusans",
    "font.size": 7,
    "axes.linewidth": 0.8,
    "axes.spines.top": False,
    "axes.spines.right": False,
})


def coordinate_edges(centers: np.ndarray) -> np.ndarray:
    """Return cell edges for a regularly spaced center coordinate."""
    centers = np.asarray(centers, dtype=float).squeeze()
    spacing = float(np.median(np.diff(centers)))
    return np.concatenate((
        [centers[0] - spacing / 2.0],
        0.5 * (centers[:-1] + centers[1:]),
        [centers[-1] + spacing / 2.0],
    ))


def interpolate_altitude(
    altitude_km: np.ndarray,
    field: np.ndarray,
    requested_altitude_km: float,
    *,
    logarithmic: bool,
) -> np.ndarray:
    """Interpolate a longitude-latitude-altitude field at one altitude."""
    altitude = np.asarray(altitude_km, dtype=float).squeeze()
    values = np.asarray(field, dtype=float)
    upper = int(np.searchsorted(altitude, requested_altitude_km))
    if upper == 0 or upper == altitude.size:
        raise ValueError("Requested altitude is outside the source grid.")
    if altitude[upper] == requested_altitude_km:
        return values[:, :, upper].copy()
    lower = upper - 1
    fraction = (
        (requested_altitude_km - altitude[lower])
        / (altitude[upper] - altitude[lower])
    )
    lower_field = values[:, :, lower]
    upper_field = values[:, :, upper]
    if logarithmic:
        lower_field = np.log(np.maximum(lower_field, 1e-300))
        upper_field = np.log(np.maximum(upper_field, 1e-300))
        return np.exp(lower_field + fraction * (upper_field - lower_field))
    return lower_field + fraction * (upper_field - lower_field)


def periodic_log_longitude_interpolation(
    source_longitude_deg: np.ndarray,
    source_field: np.ndarray,
    target_longitude_deg: np.ndarray,
) -> np.ndarray:
    """Interpolate positive longitude-latitude data across the 0/360 seam."""
    source_lon = np.asarray(source_longitude_deg, dtype=float).squeeze()
    values = np.asarray(source_field, dtype=float)
    extended_lon = np.concatenate((source_lon, [source_lon[0] + 360.0]))
    extended_log = np.concatenate((
        np.log(np.maximum(values, 1e-300)),
        np.log(np.maximum(values[:1, :], 1e-300)),
    ), axis=0)
    output = np.empty((len(target_longitude_deg), values.shape[1]))
    for latitude_index in range(values.shape[1]):
        output[:, latitude_index] = np.exp(np.interp(
            target_longitude_deg,
            extended_lon,
            extended_log[:, latitude_index],
        ))
    return output


def rounded_log_limits(field: np.ndarray) -> tuple[float, float]:
    """Return enclosing powers of ten for a positive field."""
    positive = field[np.isfinite(field) & (field > 0)]
    return (
        10.0 ** np.floor(np.log10(positive.min())),
        10.0 ** np.ceil(np.log10(positive.max())),
    )


def main() -> None:
    # F10.7 = 200 is an original source case, not the interpolated F130 case.
    gitm, amps = load_atmosphere_case(DATA, LS_DEG, F107)
    longitude = np.asarray(gitm["lon_deg"], dtype=float)
    latitude = np.asarray(gitm["lat_deg"], dtype=float)

    co2 = interpolate_altitude(
        gitm["alt_km"], gitm["nCO2"], ALTITUDE_KM, logarithmic=True
    )
    cold_o = interpolate_altitude(
        gitm["alt_km"], gitm["nO"], ALTITUDE_KM, logarithmic=True
    )
    hot_o_native = interpolate_altitude(
        amps["alt_km"], amps["nO_hot"], ALTITUDE_KM, logarithmic=True
    )
    hot_o = periodic_log_longitude_interpolation(
        amps["lon_deg"], hot_o_native, longitude
    )
    total_o = cold_o + hot_o
    temperature = interpolate_altitude(
        gitm["alt_km"], gitm["Tn"], ALTITUDE_KM, logarithmic=False
    )
    sza = interpolate_altitude(
        gitm["alt_km"], gitm["SZA_deg"], ALTITUDE_KM, logarithmic=False
    )

    longitude_edges = coordinate_edges(longitude)
    latitude_edges = coordinate_edges(latitude)
    panels = (
        (co2, r"CO$_2$ number density", "turbo",
         LogNorm(*rounded_log_limits(co2)), r"m$^{-3}$"),
        (total_o, "Total O number density", "turbo",
         LogNorm(*rounded_log_limits(total_o)), r"m$^{-3}$"),
        (temperature, "Neutral temperature", "turbo",
         mpl.colors.Normalize(temperature.min(), temperature.max()), "K"),
        (sza, "Solar zenith angle", "viridis",
         mpl.colors.Normalize(0, 180), "deg"),
    )

    fig, axes = plt.subplots(
        2, 2, figsize=(7.2, 5.3), sharex=True, sharey=True,
        constrained_layout=True,
    )
    for panel_index, (axis, panel) in enumerate(zip(axes.flat, panels)):
        field, title, cmap, norm, unit = panel
        mesh = axis.pcolormesh(
            longitude_edges, latitude_edges, field.T,
            shading="flat", cmap=cmap, norm=norm, rasterized=True,
        )
        axis.set(
            xlim=(0, 360), ylim=(-90, 90),
            xticks=(0, 90, 180, 270, 360),
            yticks=(-90, -45, 0, 45, 90),
            title=title,
        )
        axis.axvline(90, color="white", lw=0.7, ls=":", alpha=0.9)
        axis.axvline(270, color="white", lw=0.7, ls=":", alpha=0.9)
        axis.text(
            0.02, 0.96, chr(ord("a") + panel_index),
            transform=axis.transAxes, ha="left", va="top",
            fontweight="bold",
            color="white" if panel_index == 3 else "black",
        )
        axis.plot(
            2.5, -2.5, marker="*", ms=6, color="white",
            markeredgecolor="0.15", markeredgewidth=0.4,
        )
        colorbar = fig.colorbar(mesh, ax=axis, pad=0.015)
        colorbar.set_label(unit)

    axes[0, 0].set_ylabel("MSO latitude (deg)")
    axes[1, 0].set_ylabel("MSO latitude (deg)")
    axes[1, 0].set_xlabel("MSO longitude (deg)")
    axes[1, 1].set_xlabel("MSO longitude (deg)")
    fig.suptitle(
        rf"GITM + AMPS at {ALTITUDE_KM:.0f} km, "
        rf"$L_s={LS_DEG}^\circ$, F10.7 = {F107}",
        fontsize=9,
    )

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT, dpi=600, bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)
    print(f"longitude_centers_deg=({longitude[0]:g}, {longitude[-1]:g})")
    print(f"longitude_plot_limits_deg=(0, 360)")
    print(f"co2_range_m3=({co2.min():.8g}, {co2.max():.8g})")
    print(f"total_o_range_m3=({total_o.min():.8g}, {total_o.max():.8g})")
    print(
        f"temperature_range_k=({temperature.min():.8g}, "
        f"{temperature.max():.8g})"
    )
    print(f"sza_range_deg=({sza.min():.8g}, {sza.max():.8g})")
    print(f"output={OUTPUT.resolve()}")


if __name__ == "__main__":
    main()

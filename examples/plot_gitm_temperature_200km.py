"""Compare GITM neutral temperature and solar zenith angle at 200 km."""

from __future__ import annotations

from pathlib import Path
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "analysis"))
from marsaspen_analysis import load_atmosphere_case  # noqa: E402

ATMOSPHERE = REPO / "data" / "atmosphere"
OUTPUT = (
    REPO / "examples" / "figures" /
    "gitm_neutral_temperature_200km.png"
)
ALTITUDE_KM = 200.0

mpl.rcParams.update({
    "font.family": "Arial",
    "mathtext.fontset": "dejavusans",
    "font.size": 7,
    "axes.linewidth": 0.8,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "legend.frameon": False,
})


def coordinate_edges(centers: np.ndarray) -> np.ndarray:
    """Construct regularly spaced cell edges from cell centers."""
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
) -> np.ndarray:
    """Linearly interpolate a longitude-latitude-altitude field."""
    altitude = np.asarray(altitude_km, dtype=float).squeeze()
    values = np.asarray(field, dtype=float)
    if not altitude[0] <= requested_altitude_km <= altitude[-1]:
        raise ValueError(
            f"{requested_altitude_km:g} km is outside "
            f"{altitude[0]:g} to {altitude[-1]:g} km."
        )
    upper = int(np.searchsorted(altitude, requested_altitude_km))
    if altitude[upper] == requested_altitude_km:
        return values[:, :, upper].copy()
    lower = upper - 1
    fraction = (
        (requested_altitude_km - altitude[lower])
        / (altitude[upper] - altitude[lower])
    )
    return values[:, :, lower] + fraction * (
        values[:, :, upper] - values[:, :, lower]
    )


def main() -> None:
    # Use a directly supplied source case. F10.7 = 130 is not used because it
    # is interpolated between the available solar-minimum and maximum files.
    gitm, _ = load_atmosphere_case(ATMOSPHERE, 0, 200)
    temperature_lon_lat = interpolate_altitude(
        gitm["alt_km"], gitm["Tn"], ALTITUDE_KM
    )
    sza_lon_lat = interpolate_altitude(
        gitm["alt_km"], gitm["SZA_deg"], ALTITUDE_KM
    )

    # Source positions are already in MSO. Center longitude on the subsolar
    # meridian for the conventional -180 to 180 degree display.
    shift = -(temperature_lon_lat.shape[0] // 2)
    temperature_lat_lon = np.roll(temperature_lon_lat.T, shift, axis=1)
    sza_lat_lon = np.roll(sza_lon_lat.T, shift, axis=1)
    longitude_edges = np.linspace(-180.0, 180.0, 73)
    latitude_edges = coordinate_edges(gitm["lat_deg"])
    longitude_centers = np.linspace(-177.5, 177.5, 72)
    latitude_centers = np.asarray(gitm["lat_deg"], dtype=float)

    maximum_index = np.unravel_index(
        np.nanargmax(temperature_lat_lon), temperature_lat_lon.shape
    )
    maximum_longitude = longitude_centers[maximum_index[1]]
    maximum_latitude = latitude_centers[maximum_index[0]]
    minimum_sza_index = np.unravel_index(
        np.nanargmin(sza_lat_lon), sza_lat_lon.shape
    )
    minimum_sza_longitude = longitude_centers[minimum_sza_index[1]]
    minimum_sza_latitude = latitude_centers[minimum_sza_index[0]]

    fig, axes = plt.subplots(
        2, 1, figsize=(7.2, 6.2), sharex=True, sharey=True,
        constrained_layout=True,
    )
    temperature_mesh = axes[0].pcolormesh(
        longitude_edges, latitude_edges, temperature_lat_lon,
        shading="flat", cmap="turbo", rasterized=True,
    )
    sza_mesh = axes[1].pcolormesh(
        longitude_edges, latitude_edges, sza_lat_lon,
        shading="flat", cmap="viridis", vmin=0, vmax=180,
        rasterized=True,
    )

    for axis in axes:
        axis.set(
            xlim=(-180, 180), ylim=(-90, 90),
            xticks=(-180, -90, 0, 90, 180),
            yticks=(-90, -45, 0, 45, 90),
            ylabel="MSO latitude (deg)",
        )
        axis.axvline(-90, color="white", lw=0.8, ls=":", alpha=0.9)
        axis.axvline(90, color="white", lw=0.8, ls=":", alpha=0.9)

    axes[0].set_title(
        r"GITM at 200 km, $L_s=0^\circ$, F10.7 = 200", fontsize=9,
    )
    axes[1].set_xlabel("MSO longitude (deg)")
    axes[0].text(
        0.015, 0.96, "a", transform=axes[0].transAxes,
        va="top", fontweight="bold",
    )
    axes[1].text(
        0.015, 0.96, "b", transform=axes[1].transAxes,
        va="top", color="white", fontweight="bold",
    )
    axes[0].plot(
        0.0, 0.0, marker="*", ms=8, color="white",
        markeredgecolor="0.15", markeredgewidth=0.5,
        label="Subsolar point",
    )
    axes[0].plot(
        maximum_longitude, maximum_latitude,
        marker="x", ms=6, color="black", markeredgewidth=1.0,
        label=(
            "Temperature maximum\n"
            f"({maximum_longitude:.1f}°, {maximum_latitude:.1f}°)"
        ),
    )
    axes[1].plot(
        minimum_sza_longitude, minimum_sza_latitude,
        marker="*", ms=8, color="white",
        markeredgecolor="0.15", markeredgewidth=0.5,
    )
    axes[0].legend(loc="lower left", fontsize=6.5)
    temperature_colorbar = fig.colorbar(
        temperature_mesh, ax=axes[0], pad=0.015
    )
    temperature_colorbar.set_label("Neutral temperature (K)")
    sza_colorbar = fig.colorbar(sza_mesh, ax=axes[1], pad=0.015)
    sza_colorbar.set_label("Solar zenith angle (deg)")

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT, dpi=600, bbox_inches="tight")
    plt.close(fig)
    print(
        f"temperature_range_k=({temperature_lat_lon.min():.6g}, "
        f"{temperature_lat_lon.max():.6g})"
    )
    print(
        "temperature_maximum_lon_lat_deg="
        f"({maximum_longitude:.6g}, {maximum_latitude:.6g})"
    )
    print(
        "minimum_sza_lon_lat_deg="
        f"({minimum_sza_longitude:.6g}, {minimum_sza_latitude:.6g}); "
        f"sza_deg={sza_lat_lon[minimum_sza_index]:.6g}"
    )
    print(f"output={OUTPUT.resolve()}")


if __name__ == "__main__":
    main()

"""Plot MGITM neutral temperature at 200 km on its native longitude grid."""

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import loadmat

REPO = Path(__file__).resolve().parents[1]
INPUT = REPO / "data" / "atmosphere" / "gitm_ls000_f130.mat"
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


def interpolate_temperature(
    altitude_km: np.ndarray,
    temperature_k: np.ndarray,
    requested_altitude_km: float,
) -> np.ndarray:
    """Linearly interpolate the 3D neutral temperature in altitude."""
    altitude = np.asarray(altitude_km, dtype=float).squeeze()
    temperature = np.asarray(temperature_k, dtype=float)
    if not altitude[0] <= requested_altitude_km <= altitude[-1]:
        raise ValueError(
            f"{requested_altitude_km:g} km is outside "
            f"{altitude[0]:g} to {altitude[-1]:g} km."
        )
    upper = int(np.searchsorted(altitude, requested_altitude_km))
    if altitude[upper] == requested_altitude_km:
        return temperature[:, :, upper].copy()
    lower = upper - 1
    fraction = (
        (requested_altitude_km - altitude[lower])
        / (altitude[upper] - altitude[lower])
    )
    return (
        temperature[:, :, lower]
        + fraction * (
            temperature[:, :, upper] - temperature[:, :, lower]
        )
    )


def main() -> None:
    gitm = loadmat(INPUT, squeeze_me=True)
    temperature_lon_lat = interpolate_temperature(
        gitm["alt_km"], gitm["Tn"], ALTITUDE_KM
    )

    # MGITM longitude is stored from 0 to 360 degrees. Roll only for the
    # conventional -180 to 180 degree display. The packaged file does not
    # include the subsolar longitude needed to convert this native longitude
    # into MSO longitude.
    temperature_lat_lon = np.roll(
        temperature_lon_lat.T,
        -(temperature_lon_lat.shape[0] // 2),
        axis=1,
    )
    longitude_edges = np.linspace(-180.0, 180.0, 73)
    latitude_edges = coordinate_edges(gitm["lat_deg"])
    longitude_centers = np.linspace(-177.5, 177.5, 72)
    latitude_centers = np.asarray(gitm["lat_deg"], dtype=float)
    maximum_index = np.unravel_index(
        np.nanargmax(temperature_lat_lon), temperature_lat_lon.shape
    )
    maximum_longitude = longitude_centers[maximum_index[1]]
    maximum_latitude = latitude_centers[maximum_index[0]]

    fig, axis = plt.subplots(
        1, 1, figsize=(7.2, 3.5), constrained_layout=True
    )
    mesh = axis.pcolormesh(
        longitude_edges,
        latitude_edges,
        temperature_lat_lon,
        shading="flat",
        cmap="turbo",
        rasterized=True,
    )
    axis.set(
        xlim=(-180, 180),
        ylim=(-90, 90),
        xticks=(-180, -90, 0, 90, 180),
        yticks=(-90, -45, 0, 45, 90),
        xlabel="Longitude used by MarsASPEN (deg)",
        ylabel="MGITM latitude (deg)",
        title=(
            "MGITM neutral temperature at 200 km\n"
            r"$L_s=0^\circ$, F10.7 = 130"
        ),
    )
    axis.axvline(-90, color="white", lw=0.8, ls=":", alpha=0.9)
    axis.axvline(90, color="white", lw=0.8, ls=":", alpha=0.9)
    axis.plot(
        0.0, 0.0,
        marker="*", ms=8, color="white",
        markeredgecolor="0.15", markeredgewidth=0.5,
        label="MarsASPEN subsolar point",
    )
    axis.plot(
        maximum_longitude, maximum_latitude,
        marker="x", ms=6, color="black", markeredgewidth=1.0,
        label=(
            "Temperature maximum\n"
            f"({maximum_longitude:.1f}°, {maximum_latitude:.1f}°)"
        ),
    )
    axis.legend(loc="upper right", fontsize=6.5)
    colorbar = fig.colorbar(mesh, ax=axis, pad=0.015)
    colorbar.set_label("Neutral temperature (K)")

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT, dpi=600, bbox_inches="tight")
    plt.close(fig)
    print(
        f"temperature_range_k=({temperature_lat_lon.min():.6g}, "
        f"{temperature_lat_lon.max():.6g})"
    )
    print(
        f"temperature_maximum_lon_lat_deg="
        f"({maximum_longitude:.6g}, {maximum_latitude:.6g})"
    )
    print(f"output={OUTPUT.resolve()}")


if __name__ == "__main__":
    main()

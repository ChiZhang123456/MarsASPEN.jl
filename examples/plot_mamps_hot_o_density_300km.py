"""Plot the MAMPS hot-O number density at 300 km."""

from __future__ import annotations

from pathlib import Path
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "analysis"))
from marsaspen_analysis import load_atmosphere_case  # noqa: E402

ATMOSPHERE = REPO / "data" / "atmosphere"
OUTPUT = (
    REPO / "examples" / "figures" /
    "mamps_hot_o_density_300km.png"
)
ALTITUDE_KM = 300.0

mpl.rcParams.update({
    "font.family": "Arial",
    "mathtext.fontset": "dejavusans",
    "font.size": 7,
    "axes.linewidth": 0.8,
    "axes.spines.top": False,
    "axes.spines.right": False,
})


def coordinate_edges(centers: np.ndarray) -> np.ndarray:
    """Construct cell edges from regularly spaced cell centers."""
    centers = np.asarray(centers, dtype=float).squeeze()
    spacing = float(np.median(np.diff(centers)))
    return np.concatenate((
        [centers[0] - spacing / 2.0],
        0.5 * (centers[:-1] + centers[1:]),
        [centers[-1] + spacing / 2.0],
    ))


def interpolate_log_density(
    altitude_km: np.ndarray,
    density_m3: np.ndarray,
    requested_altitude_km: float,
) -> np.ndarray:
    """Log-linearly interpolate a longitude-latitude density field."""
    altitude = np.asarray(altitude_km, dtype=float).squeeze()
    density = np.asarray(density_m3, dtype=float)
    if not altitude[0] <= requested_altitude_km <= altitude[-1]:
        raise ValueError(
            f"{requested_altitude_km:g} km is outside "
            f"{altitude[0]:g} to {altitude[-1]:g} km."
        )
    upper = int(np.searchsorted(altitude, requested_altitude_km))
    if altitude[upper] == requested_altitude_km:
        return density[:, :, upper].copy()
    lower = upper - 1
    fraction = (
        (requested_altitude_km - altitude[lower])
        / (altitude[upper] - altitude[lower])
    )
    lower_log = np.log(np.maximum(density[:, :, lower], 1.0e-300))
    upper_log = np.log(np.maximum(density[:, :, upper], 1.0e-300))
    return np.exp(lower_log + fraction * (upper_log - lower_log))


def main() -> None:
    _, mamps = load_atmosphere_case(ATMOSPHERE, 0, 130)
    density_lon_lat = interpolate_log_density(
        mamps["alt_km"], mamps["nO_hot"], ALTITUDE_KM
    )

    # Preserve the native 0 to 360 degree convention. Append the periodic
    # duplicate of the 0-degree column at 360 degrees for seamless plotting.
    density_lon_lat = np.concatenate(
        (density_lon_lat, density_lon_lat[:1, :]), axis=0
    )
    longitude_centers = np.arange(0.0, 361.0, 5.0)
    density_lat_lon = density_lon_lat.T
    longitude_edges = coordinate_edges(longitude_centers)
    latitude_edges = coordinate_edges(mamps["lat_deg"])
    positive = density_lat_lon[
        np.isfinite(density_lat_lon) & (density_lat_lon > 0)
    ]
    # Rounded limits closely bracket the actual 300-km range
    # (2.49e8 to 2.30e9 m^-3) without clipping valid cells.
    limits = (2.0e8, 3.0e9)

    fig, axis = plt.subplots(
        1, 1, figsize=(7.2, 3.5), constrained_layout=True
    )
    mesh = axis.pcolormesh(
        longitude_edges,
        latitude_edges,
        np.ma.masked_less_equal(density_lat_lon, 0.0),
        shading="flat",
        cmap="turbo",
        norm=LogNorm(*limits),
        rasterized=True,
    )
    axis.set(
        xlim=(0, 360),
        ylim=(-90, 90),
        xticks=(0, 90, 180, 270, 360),
        yticks=(-90, -45, 0, 45, 90),
        xlabel="MSO longitude, native 0° to 360° convention (deg)",
        ylabel="MSO latitude (deg)",
        title=(
            "MAMPS hot-O number density at 300 km\n"
            r"$L_s=0^\circ$, F10.7 = 130"
        ),
    )
    axis.axvline(90, color="white", lw=0.8, ls=":", alpha=0.9)
    axis.axvline(270, color="white", lw=0.8, ls=":", alpha=0.9)
    axis.plot(
        0.0, 0.0, marker="*", ms=8, color="white",
        markeredgecolor="0.15", markeredgewidth=0.5,
        label="Subsolar point",
    )
    axis.legend(loc="lower left", fontsize=6.5)
    colorbar = fig.colorbar(mesh, ax=axis, pad=0.015)
    colorbar.set_ticks((2.0e8, 5.0e8, 1.0e9, 2.0e9, 3.0e9))
    colorbar.set_ticklabels(
        (r"$2\times10^8$", r"$5\times10^8$", r"$10^9$",
         r"$2\times10^9$", r"$3\times10^9$")
    )
    colorbar.set_label(r"Hot-O number density (m$^{-3}$)")

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT, dpi=600, bbox_inches="tight")
    plt.close(fig)
    print(
        f"density_range_m3=({positive.min():.9g}, "
        f"{positive.max():.9g})"
    )
    print(f"color_limits_m3={limits}")
    print(f"output={OUTPUT.resolve()}")


if __name__ == "__main__":
    main()

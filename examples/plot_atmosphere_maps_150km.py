"""Plot global CO2 and total-O density maps at 150 km.

The figure contains every packaged MGITM and MAMPS atmosphere case. Rows are
organized first by species and then by season, while columns are the three
F10.7 values. CO2 is taken from MGITM. Total O is the sum of MGITM cold O and
MAMPS hot O after each field is interpolated logarithmically to 150 km.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
from scipy.io import loadmat

REPO = Path(__file__).resolve().parents[1]
DATA = REPO / "data" / "atmosphere"
OUTPUT = (
    REPO / "examples" / "figures" /
    "gitm_mamps_co2_o_maps_150km.png"
)

LS_VALUES = (0, 90, 180, 270)
F107_VALUES = (70, 130, 200)
ALTITUDE_KM = 150.0

mpl.rcParams.update({
    "font.family": "Arial",
    "mathtext.fontset": "dejavusans",
    "font.size": 7,
    "axes.linewidth": 0.8,
    "axes.spines.top": False,
    "axes.spines.right": False,
})


def interpolate_log_altitude(
    altitude_km: np.ndarray,
    density_m3: np.ndarray,
    requested_altitude_km: float,
) -> np.ndarray:
    """Interpolate a longitude-latitude density field in log density."""
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


def coordinate_edges(centers: np.ndarray, periodic: bool) -> np.ndarray:
    """Construct cell edges from regularly spaced coordinate centers."""
    centers = np.asarray(centers, dtype=float).squeeze()
    spacing = float(np.median(np.diff(centers)))
    edges = np.concatenate((
        [centers[0] - spacing / 2],
        0.5 * (centers[:-1] + centers[1:]),
        [centers[-1] + spacing / 2],
    ))
    if periodic:
        edges[0], edges[-1] = 0.0, 360.0
    return edges


def rounded_log_limits(arrays: list[np.ndarray]) -> tuple[float, float]:
    """Return common decade limits covering every positive array value."""
    positive = np.concatenate([
        values[np.isfinite(values) & (values > 0)]
        for values in arrays
    ])
    return (
        10.0 ** np.floor(np.log10(positive.min())),
        10.0 ** np.ceil(np.log10(positive.max())),
    )


def main() -> None:
    cases: dict[tuple[int, int], dict[str, np.ndarray]] = {}
    co2_fields: list[np.ndarray] = []
    oxygen_fields: list[np.ndarray] = []

    for ls in LS_VALUES:
        for f107 in F107_VALUES:
            suffix = f"ls{ls:03d}_f{f107:03d}.mat"
            gitm = loadmat(DATA / f"gitm_{suffix}", squeeze_me=True)
            mamps = loadmat(DATA / f"mamps_{suffix}", squeeze_me=True)

            co2 = interpolate_log_altitude(
                gitm["alt_km"], gitm["nCO2"], ALTITUDE_KM
            )
            cold_o = interpolate_log_altitude(
                gitm["alt_km"], gitm["nO"], ALTITUDE_KM
            )
            hot_o = interpolate_log_altitude(
                mamps["alt_km"], mamps["nO_hot"], ALTITUDE_KM
            )
            total_o = cold_o + hot_o
            cases[(ls, f107)] = {
                "longitude": np.asarray(gitm["lon_deg"], dtype=float),
                "latitude": np.asarray(gitm["lat_deg"], dtype=float),
                "co2": co2,
                "oxygen": total_o,
            }
            co2_fields.append(co2)
            oxygen_fields.append(total_o)

    co2_limits = rounded_log_limits(co2_fields)
    oxygen_limits = rounded_log_limits(oxygen_fields)
    norms = {
        "co2": LogNorm(*co2_limits),
        "oxygen": LogNorm(*oxygen_limits),
    }

    fig, axes = plt.subplots(
        8, 3, figsize=(7.2, 11.5),
        sharex=True, sharey=True, constrained_layout=True,
    )
    last_mesh: dict[str, mpl.collections.Collection] = {}
    panel_index = 0

    for species_block, (field, species_label) in enumerate((
        ("co2", r"CO$_2$"),
        ("oxygen", "Total O"),
    )):
        for season_row, ls in enumerate(LS_VALUES):
            row = species_block * len(LS_VALUES) + season_row
            for col, f107 in enumerate(F107_VALUES):
                axis = axes[row, col]
                case = cases[(ls, f107)]
                longitude_edges = coordinate_edges(
                    case["longitude"], periodic=True
                )
                latitude_edges = coordinate_edges(
                    case["latitude"], periodic=False
                )
                mesh = axis.pcolormesh(
                    longitude_edges,
                    latitude_edges,
                    case[field].T,
                    cmap="turbo",
                    norm=norms[field],
                    shading="flat",
                    rasterized=True,
                )
                last_mesh[field] = mesh
                axis.set_xlim(0, 360)
                axis.set_ylim(-90, 90)
                axis.set_xticks((0, 90, 180, 270, 360))
                axis.set_yticks((-90, -45, 0, 45, 90))
                axis.text(
                    0.02, 0.96, chr(ord("a") + panel_index),
                    transform=axis.transAxes, ha="left", va="top",
                    fontsize=7.5, fontweight="bold", color="black",
                )
                panel_index += 1
                if row == 0:
                    axis.set_title(f"F10.7 = {f107}", fontsize=8)
                if col == 0:
                    axis.set_ylabel(
                        species_label + "\n"
                        + rf"$L_s={ls}^\circ$"
                        + "\nLatitude (deg)"
                    )
                if row == 7:
                    axis.set_xlabel("Longitude (deg)")

    co2_colorbar = fig.colorbar(
        last_mesh["co2"], ax=axes[:4, :],
        location="right", pad=0.015, shrink=0.82,
    )
    co2_colorbar.set_label(
        r"CO$_2$ number density (m$^{-3}$)"
    )
    oxygen_colorbar = fig.colorbar(
        last_mesh["oxygen"], ax=axes[4:, :],
        location="right", pad=0.015, shrink=0.82,
    )
    oxygen_colorbar.set_label(
        r"Total O number density (m$^{-3}$)"
    )
    fig.suptitle(
        "MGITM CO$_2$ and MGITM + MAMPS total O at 150 km",
        fontsize=9,
    )
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT, dpi=600, bbox_inches="tight", pad_inches=0.08)
    plt.close(fig)
    print(f"CO2_limits_m3={co2_limits}")
    print(f"total_O_limits_m3={oxygen_limits}")
    print(f"output={OUTPUT.resolve()}")


if __name__ == "__main__":
    main()

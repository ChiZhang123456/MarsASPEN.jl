"""Plot all packaged MGITM and MAMPS density profiles.

The 4 by 3 panel grid spans four seasons and three solar-activity levels.
Every profile is sampled at 0 degrees longitude and 0 degrees latitude.
MGITM contains the cold neutral species CO2, O, O2, N2, and CO. MAMPS supplies
the separate hot-O population. Number densities are displayed in SI units.

MarsASPEN extrapolates the first two logarithmic density levels down to 80 km.
Above the native MGITM top boundary, each cold species is hydrostatically
extrapolated using the local top-layer neutral temperature. MAMPS hot O is
used only inside its native altitude range.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import loadmat

REPO = Path(__file__).resolve().parents[1]
DATA = REPO / "data" / "atmosphere"
OUTPUT = REPO / "examples" / "figures" / "gitm_mamps_density_cases.png"

LS_VALUES = (0, 90, 180, 270)
F107_VALUES = (70, 130, 200)
SPECIES = (
    ("nCO2", r"CO$_2$", "#4C72B0", 44.01),
    ("nO", "O", "#55A868", 15.999),
    ("nO2", r"O$_2$", "#C44E52", 31.998),
    ("nN2", r"N$_2$", "#8172B2", 28.014),
    ("nCO", "CO", "#937860", 28.010),
)
AMU_KG = 1.66053906660e-27
KB_J_K = 1.380649e-23
MARS_G0_M_S2 = 3.71
MARS_RADIUS_KM = 3388.25

mpl.rcParams.update({
    "font.family": "Arial",
    "mathtext.fontset": "dejavusans",
    "font.size": 7,
    "axes.linewidth": 0.8,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "legend.frameon": False,
})


def horizontal_profile(
    data: dict[str, np.ndarray], field: str, logarithmic: bool,
) -> np.ndarray:
    """Bilinearly interpolate one vertical column at lon=0, lat=0."""
    lon = np.asarray(data["lon_deg"]).squeeze()
    lat = np.asarray(data["lat_deg"]).squeeze()
    # Longitude zero lies across the periodic boundary between the last and
    # first grid cells. Latitude zero lies between -2.5 and +2.5 degrees.
    i0, i1 = lon.size - 1, 0
    wx = (360.0 - lon[i0]) / (lon[i1] + 360.0 - lon[i0])
    j0 = int(np.searchsorted(lat, 0.0) - 1)
    j1 = j0 + 1
    wy = (0.0 - lat[j0]) / (lat[j1] - lat[j0])
    values = np.asarray(data[field], dtype=float)
    if logarithmic:
        values = np.log(np.maximum(values, 1.0e-300))
    c0 = (1.0 - wx) * values[i0, j0, :] + wx * values[i1, j0, :]
    c1 = (1.0 - wx) * values[i0, j1, :] + wx * values[i1, j1, :]
    profile = (1.0 - wy) * c0 + wy * c1
    return np.exp(profile) if logarithmic else profile


def cold_model_profile(
    altitude_native: np.ndarray,
    density_native: np.ndarray,
    top_temperature_k: float,
    mass_amu: float,
    altitude_plot: np.ndarray,
) -> np.ndarray:
    """Match the model's lower and hydrostatic upper extrapolations."""
    log_density = np.log(np.maximum(density_native, 1.0e-300))
    result = np.empty(altitude_plot.shape)
    inside = altitude_plot <= altitude_native[-1]
    result[inside] = np.exp(np.interp(
        altitude_plot[inside], altitude_native, log_density
    ))
    below = altitude_plot < altitude_native[0]
    slope = (
        (log_density[1] - log_density[0])
        / (altitude_native[1] - altitude_native[0])
    )
    result[below] = np.exp(
        log_density[0] + slope * (altitude_plot[below] - altitude_native[0])
    )
    above = altitude_plot > altitude_native[-1]
    gravity = MARS_G0_M_S2 * (
        MARS_RADIUS_KM / (MARS_RADIUS_KM + altitude_native[-1])
    ) ** 2
    scale_height_km = (
        KB_J_K * top_temperature_k / (mass_amu * AMU_KG * gravity) / 1000.0
    )
    result[above] = density_native[-1] * np.exp(
        -(altitude_plot[above] - altitude_native[-1]) / scale_height_km
    )
    return result


def hot_model_profile(
    altitude_native: np.ndarray,
    density_native: np.ndarray,
    altitude_plot: np.ndarray,
) -> np.ndarray:
    """Interpolate MAMPS hot O and mask points outside its native range."""
    result = np.full(altitude_plot.shape, np.nan)
    inside = (
        (altitude_plot >= altitude_native[0])
        & (altitude_plot <= altitude_native[-1])
    )
    result[inside] = np.exp(np.interp(
        altitude_plot[inside], altitude_native,
        np.log(np.maximum(density_native, 1.0e-300)),
    ))
    return result


def main() -> None:
    altitude = np.arange(80.0, 1001.0, 1.0)
    fig, axes = plt.subplots(
        len(LS_VALUES), len(F107_VALUES), figsize=(7.2, 8.5),
        sharex=True, sharey=True, constrained_layout=True,
    )

    for row, ls in enumerate(LS_VALUES):
        for col, f107 in enumerate(F107_VALUES):
            axis = axes[row, col]
            suffix = f"ls{ls:03d}_f{f107:03d}.mat"
            gitm = loadmat(DATA / f"gitm_{suffix}", squeeze_me=True)
            mamps = loadmat(DATA / f"mamps_{suffix}", squeeze_me=True)
            gitm_alt = np.asarray(gitm["alt_km"], dtype=float).squeeze()

            # Plot each cold neutral only over the vertical range supported by
            # MGITM, plus the model's prescribed lower extrapolation to 80 km.
            temperature_profile = horizontal_profile(gitm, "Tn", False)
            top_temperature = float(temperature_profile[-1])
            for field, label, color, mass_amu in SPECIES:
                native = horizontal_profile(gitm, field, True)
                profile = cold_model_profile(
                    gitm_alt, native, top_temperature, mass_amu, altitude
                )
                axis.plot(profile, altitude, color=color, lw=1.0, label=label)

            # MAMPS uses its own altitude grid. Hot O is visually separated by
            # a dashed dark-green curve.
            hot_alt = np.asarray(mamps["alt_km"], dtype=float).squeeze()
            hot_native = horizontal_profile(mamps, "nO_hot", True)
            hot_profile = hot_model_profile(hot_alt, hot_native, altitude)
            axis.plot(
                hot_profile, altitude, color="#006D5B", lw=1.2, ls="--",
                label="O hot (MAMPS)",
            )

            axis.set_xscale("log")
            axis.set_xlim(1.0e0, 1.0e19)
            axis.set_ylim(80, 1000)
            axis.grid(True, which="major", color="0.90", lw=0.5)
            axis.set_title(rf"$L_s={ls}^\circ$, F10.7={f107}")
            axis.text(
                0.02, 0.97, chr(ord("a") + row * 3 + col),
                transform=axis.transAxes, ha="left", va="top",
                fontweight="bold", fontsize=8,
            )

    for axis in axes[-1, :]:
        axis.set_xlabel(r"Number density (m$^{-3}$)")
    for axis in axes[:, 0]:
        axis.set_ylabel("Altitude (km)")
    axes[0, 0].legend(
        ncol=2, fontsize=6, loc="lower left", bbox_to_anchor=(0.0, 1.18)
    )
    fig.suptitle(
        r"MGITM cold neutrals and MAMPS hot O at lon=0$^\circ$, lat=0$^\circ$",
        fontsize=9,
    )
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT, dpi=600, bbox_inches="tight")
    plt.close(fig)
    print(f"output={OUTPUT.resolve()}")


if __name__ == "__main__":
    main()

"""Plot H+ and H-ENA density profiles from 100 to 600 km.

The profiles are global spherical-area means calculated from the 5 degree by
5 degree by 1 km spatial moments produced by the ten-million-particle uniform
dayside H+ simulation. The second panel shows the neutral fraction, which
quantifies where charge-exchanged H-ENA becomes the dominant projectile state.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "analysis"))

from marsaspen_analysis import load_history_mat  # noqa: E402

mpl.rcParams.update({
    "font.family": "Arial",
    "mathtext.fontset": "dejavusans",
    "font.size": 7,
    "axes.linewidth": 0.8,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "legend.frameon": False,
})

ALTITUDE_MIN_KM = 100.0
ALTITUDE_MAX_KM = 600.0


def spherical_area_mean(
    values: np.ndarray,
    latitude_edges_deg: np.ndarray,
    longitude_centers_deg: np.ndarray,
) -> np.ndarray:
    """Return the spherical-area-weighted horizontal mean at each altitude."""
    latitude_weights = np.diff(np.sin(np.deg2rad(latitude_edges_deg)))
    weights = latitude_weights[:, None] * np.ones(
        (1, longitude_centers_deg.size)
    )
    return (
        np.sum(values * weights[None, :, :], axis=(1, 2)) /
        np.sum(weights)
    )


def interpolated_crossing(
    altitude_km: np.ndarray,
    hplus_density: np.ndarray,
    hena_density: np.ndarray,
) -> float:
    """Find the altitude where the two positive density profiles are equal."""
    valid = (
        (hplus_density > 0.0) &
        (hena_density > 0.0) &
        np.isfinite(hplus_density) &
        np.isfinite(hena_density)
    )
    log_ratio = np.full_like(hplus_density, np.nan)
    log_ratio[valid] = np.log10(hena_density[valid] / hplus_density[valid])
    crossing_indices = np.flatnonzero(
        np.isfinite(log_ratio[:-1]) &
        np.isfinite(log_ratio[1:]) &
        (log_ratio[:-1] * log_ratio[1:] <= 0.0)
    )
    if crossing_indices.size == 0:
        return float("nan")
    index = int(crossing_indices[0])
    fraction = (
        -log_ratio[index] /
        (log_ratio[index + 1] - log_ratio[index])
    )
    return float(
        altitude_km[index] +
        fraction * (altitude_km[index + 1] - altitude_km[index])
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "moments_file", nargs="?", type=Path,
        default=REPO / "examples" / "output" /
        "dayside_hplus_10000000_3d_moments.mat",
    )
    parser.add_argument(
        "--output", type=Path,
        default=REPO / "examples" / "figures" /
        "dayside_hplus_10000000_density_100_600km.png",
    )
    args = parser.parse_args()

    moments = load_history_mat(args.moments_file)
    altitude = np.ravel(moments["altitude_centers_km"])
    latitude_edges = np.ravel(moments["latitude_edges_deg"])
    longitude_centers = np.ravel(moments["longitude_centers_deg"])
    density = np.asarray(
        moments["number_density_by_charge_m3"], dtype=float
    )
    downward_flux = np.asarray(
        moments["downward_radial_flux_by_charge_m2_s"], dtype=float
    )
    upward_flux = np.asarray(
        moments["upward_radial_flux_by_charge_m2_s"], dtype=float
    )

    mean = lambda values: spherical_area_mean(  # noqa: E731
        values, latitude_edges, longitude_centers
    )
    hena_density = mean(density[0])
    hplus_density = mean(density[1])
    hena_downward_flux = mean(downward_flux[0])
    hena_upward_flux = mean(upward_flux[0])
    total_density = hplus_density + hena_density
    neutral_fraction = np.divide(
        hena_density,
        total_density,
        out=np.full_like(total_density, np.nan),
        where=total_density > 0.0,
    )
    crossing_altitude = interpolated_crossing(
        altitude, hplus_density, hena_density
    )

    selected = (
        (altitude >= ALTITUDE_MIN_KM) &
        (altitude < ALTITUDE_MAX_KM)
    )
    plotted_altitude = altitude[selected]

    fig, axes = plt.subplots(
        1, 2, figsize=(7.2, 4.5), sharey=True, layout="constrained"
    )
    axes[0].plot(
        np.ma.masked_less_equal(hplus_density[selected], 0.0),
        plotted_altitude,
        color="#C44E52", lw=1.5, label="H$^+$",
    )
    axes[0].plot(
        np.ma.masked_less_equal(hena_density[selected], 0.0),
        plotted_altitude,
        color="#4C72B0", lw=1.5, label="H-ENA",
    )
    axes[0].set(
        xscale="log",
        xlabel=r"Number density (m$^{-3}$)",
        ylabel="Altitude (km)",
        title="Global spherical-mean density",
        ylim=(ALTITUDE_MIN_KM, ALTITUDE_MAX_KM),
    )
    axes[0].legend(loc="lower right")

    axes[1].plot(
        neutral_fraction[selected],
        plotted_altitude,
        color="#4C72B0", lw=1.5,
    )
    axes[1].axvline(0.5, color="0.35", lw=0.8, ls=":")
    axes[1].set(
        xlim=(0.0, 1.0),
        xlabel=r"Neutral fraction, $n_{\mathrm{H\!-\!ENA}}/"
        r"(n_{\mathrm{H^+}}+n_{\mathrm{H\!-\!ENA}})$",
        title="H-ENA fraction",
    )

    if np.isfinite(crossing_altitude):
        for axis in axes:
            axis.axhline(
                crossing_altitude, color="0.35", lw=0.8, ls="--"
            )
        axes[1].annotate(
            f"Equal densities\n{crossing_altitude:.1f} km",
            xy=(0.5, crossing_altitude),
            xytext=(0.62, crossing_altitude + 35.0),
            arrowprops={"arrowstyle": "-", "color": "0.3", "lw": 0.7},
            ha="left", va="center",
        )

    for panel, axis in enumerate(axes):
        axis.grid(True, which="major", color="0.90", lw=0.5)
        axis.text(
            0.025, 0.975, "ab"[panel],
            transform=axis.transAxes, ha="left", va="top",
            fontsize=8, fontweight="bold",
        )

    fig.suptitle(
        "H$^+$ to H-ENA conversion in the uniform-dayside simulation\n"
        r"10,000,000 injected H$^+$, $\mathbf{U}=(-400,0,0)$ km/s, "
        r"$T=10$ eV, $n=5$ cm$^{-3}$",
        fontsize=9,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=600, bbox_inches="tight")
    plt.close(fig)

    print(f"output={args.output.resolve()}")
    print(f"equal_density_altitude_km={crossing_altitude:.6g}")
    for target_altitude in (120.5, 150.5, 200.5, 250.5, 300.5, 400.5,
                            500.5, 599.5):
        index = int(np.argmin(np.abs(altitude - target_altitude)))
        print(
            f"altitude_km={altitude[index]:.1f} "
            f"hplus_density_m3={hplus_density[index]:.9g} "
            f"hena_density_m3={hena_density[index]:.9g} "
            f"hena_fraction={neutral_fraction[index]:.6f}"
        )
    top_index = int(np.argmin(np.abs(altitude - 599.5)))
    print(
        f"at_599p5km_hena_downward_flux_m2_s="
        f"{hena_downward_flux[top_index]:.9g}"
    )
    print(
        f"at_599p5km_hena_upward_flux_m2_s="
        f"{hena_upward_flux[top_index]:.9g}"
    )


if __name__ == "__main__":
    main()

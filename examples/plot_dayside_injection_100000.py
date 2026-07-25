"""Plot the position and MSO velocity distributions of a dayside injector."""

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


def gaussian(x: np.ndarray, mean: float, sigma: float) -> np.ndarray:
    """Normalized one-dimensional Gaussian probability density."""
    return np.exp(-0.5 * ((x - mean) / sigma) ** 2) / (
        np.sqrt(2.0 * np.pi) * sigma
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input_mat", type=Path)
    parser.add_argument(
        "--output", type=Path,
        default=REPO / "examples" / "figures" /
        "dayside_hplus_injection_100000_4panel.png",
    )
    args = parser.parse_args()
    data = load_history_mat(args.input_mat)
    longitude = np.ravel(data["longitude_deg"])
    latitude = np.ravel(data["latitude_deg"])
    velocity = np.asarray(data["velocity_m_s"], dtype=float) / 1000.0
    if velocity.shape[0] == 3 and velocity.shape[1] == longitude.size:
        velocity = velocity.T
    weights = np.ravel(data["density_weight_m3"])

    fig, axes = plt.subplots(
        2, 2, figsize=(7.2, 5.4), constrained_layout=True
    )

    # Divide raw bin counts by spherical-bin solid angle. A uniform surface
    # distribution should then appear spatially uniform despite convergence
    # of longitude lines toward the poles.
    lon_edges = np.linspace(-90.0, 90.0, 73)
    lat_edges = np.linspace(-90.0, 90.0, 37)
    counts, _, _ = np.histogram2d(
        longitude, latitude, bins=(lon_edges, lat_edges)
    )
    dlon = np.deg2rad(np.diff(lon_edges))[:, None]
    dsinlat = np.diff(np.sin(np.deg2rad(lat_edges)))[None, :]
    particles_per_sr = counts / (dlon * dsinlat)
    mesh = axes[0, 0].pcolormesh(
        lon_edges, lat_edges, particles_per_sr.T,
        cmap="turbo", shading="flat", rasterized=True,
    )
    axes[0, 0].set(
        xlim=(-90, 90), ylim=(-90, 90),
        xlabel="MSO longitude (deg)", ylabel="MSO latitude (deg)",
        title="Initial position at 600 km",
    )
    colorbar = fig.colorbar(mesh, ax=axes[0, 0], pad=0.02)
    colorbar.set_label(r"Sample density (particles sr$^{-1}$)")

    labels = (r"$V_x$", r"$V_y$", r"$V_z$")
    colors = ("#C44E52", "#4C72B0", "#55A868")
    expected_means = (-400.0, 0.0, 0.0)
    mass_hplus = 1.007276466621 * 1.66053906660e-27
    sigma_expected = np.sqrt(
        10.0 * 1.602176634e-19 / mass_hplus
    ) / 1000.0
    velocity_axes = (axes[0, 1], axes[1, 0], axes[1, 1])
    for component, axis in enumerate(velocity_axes):
        values = velocity[:, component]
        mean = np.average(values, weights=weights)
        sigma = np.sqrt(np.average((values - mean) ** 2, weights=weights))
        center = expected_means[component]
        bounds = (center - 5 * sigma_expected, center + 5 * sigma_expected)
        bins = np.linspace(*bounds, 61)
        axis.hist(
            values, bins=bins, weights=weights, density=True,
            histtype="stepfilled", alpha=0.68, color=colors[component],
            edgecolor=colors[component], linewidth=0.8,
        )
        x = np.linspace(*bounds, 400)
        axis.plot(
            x, gaussian(x, center, sigma_expected),
            color="0.15", lw=1.0, ls="--", label="10 eV Maxwellian",
        )
        axis.set(
            xlim=bounds,
            xlabel=f"{labels[component]} (km/s)",
            ylabel="Probability density",
            title=f"{labels[component]} distribution",
        )
        axis.text(
            0.97, 0.95,
            f"mean = {mean:.2f} km/s\nstd = {sigma:.2f} km/s",
            transform=axis.transAxes, ha="right", va="top",
        )
        axis.legend(loc="upper left", fontsize=6.5)

    for label, axis in zip("abcd", axes.flat):
        axis.text(
            0.025, 0.975, label, transform=axis.transAxes,
            ha="left", va="top", fontsize=8, fontweight="bold",
        )
    fig.suptitle(
        "Uniform dayside H$^+$ injection\n"
        r"100,000 macro particles, 600 km, "
        r"$\mathbf{U}=(-400,0,0)$ km/s, $T=10$ eV, "
        r"$n=5$ cm$^{-3}$",
        fontsize=9,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=600, bbox_inches="tight")
    plt.close(fig)
    print(f"output={args.output.resolve()}")


if __name__ == "__main__":
    main()

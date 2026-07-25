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
    position = np.asarray(data["position_m"], dtype=float) / 1000.0
    if position.shape[0] == 3 and position.shape[1] == longitude.size:
        position = position.T
    velocity = np.asarray(data["velocity_m_s"], dtype=float) / 1000.0
    if velocity.shape[0] == 3 and velocity.shape[1] == longitude.size:
        velocity = velocity.T
    weights = np.ravel(data["density_weight_m3"])

    fig = plt.figure(figsize=(7.2, 5.4), constrained_layout=True)
    position_axis = fig.add_subplot(2, 2, 1, projection="3d")
    velocity_axes = (
        fig.add_subplot(2, 2, 2),
        fig.add_subplot(2, 2, 3),
        fig.add_subplot(2, 2, 4),
    )

    # Draw Mars at the origin and a deterministic subset of the 100,000
    # injection positions on the 600 km dayside spherical surface.
    mars_radius_km = 3388.25
    theta = np.linspace(0.0, 2.0 * np.pi, 80)
    polar = np.linspace(0.0, np.pi, 40)
    sphere_x = mars_radius_km * np.outer(
        np.cos(theta), np.sin(polar)
    )
    sphere_y = mars_radius_km * np.outer(
        np.sin(theta), np.sin(polar)
    )
    sphere_z = mars_radius_km * np.outer(
        np.ones_like(theta), np.cos(polar)
    )
    position_axis.plot_surface(
        sphere_x, sphere_y, sphere_z, color="#C96D3B",
        alpha=0.82, linewidth=0, shade=True, rasterized=True,
    )
    rng = np.random.default_rng(61)
    selected = rng.choice(position.shape[0], size=3000, replace=False)
    position_axis.scatter(
        position[selected, 0], position[selected, 1], position[selected, 2],
        s=1.0, color="#277DA1", alpha=0.28, depthshade=False,
        rasterized=True, label="H$^+$ at 600 km",
    )
    limit = mars_radius_km + 900.0
    position_axis.set(
        xlim=(-limit, limit), ylim=(-limit, limit), zlim=(-limit, limit),
        xlabel="MSO X (km)", ylabel="MSO Y (km)", zlabel="MSO Z (km)",
        title="Initial positions on the dayside",
    )
    position_axis.set_box_aspect((1, 1, 1))
    # Look approximately from the Sun toward Mars so the selected dayside
    # hemisphere and its 600 km stand-off from the planet are immediately
    # visible.
    position_axis.view_init(elev=18, azim=0)
    position_axis.legend(loc="upper left", fontsize=6.5)
    position_axis.text2D(
        0.69, 0.91, "Sunward\n+X", transform=position_axis.transAxes,
        ha="center", va="center", color="#B8860B", fontweight="bold",
    )
    position_axis.quiver(
        0, 0, 0, limit, 0, 0, color="#E3B341",
        arrow_length_ratio=0.08, linewidth=1.2,
    )

    labels = (r"$V_x$", r"$V_y$", r"$V_z$")
    colors = ("#C44E52", "#4C72B0", "#55A868")
    expected_means = (-400.0, 0.0, 0.0)
    mass_hplus = 1.007276466621 * 1.66053906660e-27
    sigma_expected = np.sqrt(
        10.0 * 1.602176634e-19 / mass_hplus
    ) / 1000.0
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

    position_axis.text2D(
        0.025, 0.975, "a", transform=position_axis.transAxes,
        ha="left", va="top", fontsize=8, fontweight="bold",
    )
    for label, axis in zip("bcd", velocity_axes):
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

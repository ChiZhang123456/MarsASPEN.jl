"""Plot the MarsASPEN scattering-angle inverse CDF and random samples."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

REPO = Path(__file__).resolve().parents[1]

mpl.rcParams.update(
    {
        "font.family": "Arial",
        "mathtext.fontset": "dejavusans",
        "font.size": 8,
        "axes.linewidth": 0.8,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "legend.frameon": False,
        "xtick.direction": "out",
        "ytick.direction": "out",
    }
)


def load_inverse_cdf(filename: Path) -> tuple[np.ndarray, np.ndarray]:
    """Read random-number and laboratory-angle columns from the table."""
    # The packaged table begins with six comment lines followed by one
    # human-readable, non-comment column-header line.
    table = np.loadtxt(filename, comments="#", skiprows=7, dtype=float)
    probability = np.asarray(table[:, 0], dtype=float)
    angle_deg = np.asarray(table[:, 1], dtype=float)
    order = np.argsort(probability)
    return probability[order], angle_deg[order]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", type=int, default=1_000_000)
    parser.add_argument("--seed", type=int, default=73)
    parser.add_argument(
        "--output",
        type=Path,
        default=REPO / "examples" / "figures" /
        "scattering_angle_sampling.png",
    )
    args = parser.parse_args()

    inverse_cdf_file = (
        REPO / "data" / "cross_sections" /
        "scattering_angle_distribution.txt"
    )
    probability, angle_deg = load_inverse_cdf(inverse_cdf_file)

    # MarsASPEN uses the same inverse-transform operation in Julia:
    # theta = interp(U_theta, random_grid, angle_grid).
    rng = np.random.default_rng(args.seed)
    uniform_random = rng.random(args.samples)
    sampled_angle = np.interp(uniform_random, probability, angle_deg)

    # Logarithmic angle bins resolve the dominant small-angle population and
    # the rare large-angle tail in the same panel.
    bin_edges = np.geomspace(max(angle_deg[0], 0.1), angle_deg[-1], 61)
    sampled_count, _ = np.histogram(sampled_angle, bins=bin_edges)
    sampled_fraction = sampled_count / args.samples
    expected_cdf = np.interp(bin_edges, angle_deg, probability)
    expected_fraction = np.diff(expected_cdf)
    bin_center = np.sqrt(bin_edges[:-1] * bin_edges[1:])

    figure, axes = plt.subplots(
        1, 2, figsize=(7.2, 3.15), constrained_layout=True
    )

    axes[0].plot(
        probability, angle_deg, color="#2878B5", linewidth=1.3
    )
    axes[0].set(
        xlim=(0, 1),
        yscale="log",
        ylim=(0.1, 200),
        xlabel=r"Uniform random number, $U_\theta$",
        ylabel=r"Laboratory scattering angle, $\theta$ (deg)",
        title="Inverse-CDF lookup table",
    )

    axes[1].stairs(
        sampled_fraction,
        bin_edges,
        fill=True,
        color="#E07B39",
        alpha=0.60,
        linewidth=0.8,
        label=f"Random samples, N = {args.samples:,}",
    )
    axes[1].plot(
        bin_center,
        expected_fraction,
        color="#222222",
        linewidth=1.0,
        linestyle="--",
        label="Probability from lookup table",
    )
    positive = np.r_[sampled_fraction[sampled_fraction > 0],
                     expected_fraction[expected_fraction > 0]]
    lower_limit = max(np.min(positive) * 0.5, 1e-8)
    axes[1].set(
        xscale="log",
        yscale="log",
        xlim=(0.1, 200),
        ylim=(lower_limit, min(1.0, np.max(positive) * 2.0)),
        xlabel=r"Laboratory scattering angle, $\theta$ (deg)",
        ylabel="Probability per logarithmic bin",
        title="Monte Carlo scattering-angle distribution",
    )
    axes[1].legend(loc="best", fontsize=7)
    axes[1].text(
        0.97,
        0.05,
        "Median = "
        f"{np.median(sampled_angle):.2f} deg\n"
        "95th percentile = "
        f"{np.quantile(sampled_angle, 0.95):.2f} deg\n"
        "99th percentile = "
        f"{np.quantile(sampled_angle, 0.99):.1f} deg",
        transform=axes[1].transAxes,
        ha="right",
        va="bottom",
        fontsize=7,
        bbox={
            "boxstyle": "round,pad=0.25",
            "facecolor": "white",
            "edgecolor": "0.75",
            "alpha": 0.92,
        },
    )

    for label, axis in zip(("a", "b"), axes):
        axis.text(
            0.02,
            0.97,
            label,
            transform=axis.transAxes,
            ha="left",
            va="top",
            fontsize=9,
            fontweight="bold",
        )
        axis.grid(True, which="both", color="0.91", linewidth=0.5)

    figure.suptitle(
        "Random scattering-angle sampling used by MarsASPEN",
        fontsize=10,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.output, dpi=600, bbox_inches="tight")
    plt.close(figure)

    print(f"median_angle_deg={np.median(sampled_angle):.6g}")
    print(f"p95_angle_deg={np.quantile(sampled_angle, 0.95):.6g}")
    print(f"p99_angle_deg={np.quantile(sampled_angle, 0.99):.6g}")
    print(f"output={args.output.resolve()}")


if __name__ == "__main__":
    main()

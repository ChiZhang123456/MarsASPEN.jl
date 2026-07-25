"""Plot O and CO2 ionization rates for H ENA and H+ source simulations.

The 2 by 2 quantitative grid uses target species as rows and source species as
columns. Every panel shows the instantaneous H ENA projectile contribution,
the H+ projectile contribution, and their sum from 100 to 400 km.
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

from marsaspen_analysis import ionization_rate_from_mat, load_history_mat  # noqa: E402

mpl.rcParams.update({
    "font.family": "Arial",
    "mathtext.fontset": "dejavusans",
    "font.size": 7,
    "axes.linewidth": 0.8,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "legend.frameon": False,
})

COLORS = ("#4C72B0", "#C44E52")
CHARGE_LABELS = ("H ENA contribution", r"H$^+$ contribution")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("oxygen_h_ena", type=Path)
    parser.add_argument("oxygen_hplus", type=Path)
    parser.add_argument("co2_h_ena", type=Path)
    parser.add_argument("co2_hplus", type=Path)
    parser.add_argument(
        "--output",
        type=Path,
        default=REPO / "examples" / "figures" /
        "oxygen_co2_ionization_rate_4panel.png",
    )
    args = parser.parse_args()

    paths = (
        (args.oxygen_h_ena, args.oxygen_hplus),
        (args.co2_h_ena, args.co2_hplus),
    )
    results = [
        [ionization_rate_from_mat(load_history_mat(path)) for path in row]
        for row in paths
    ]

    fig, axes = plt.subplots(
        2, 2, figsize=(7.2, 5.5), sharey=True, constrained_layout=True
    )
    target_labels = ("O", r"CO$_2$")
    source_titles = ("H ENA source", r"H$^+$ source")
    panel_labels = ("a", "b", "c", "d")

    for row in range(2):
        for col in range(2):
            axis = axes[row, col]
            result = results[row][col]
            altitude = result["altitude_km"]
            by_charge = result["rate_by_charge_m3_s1"]
            total = result["total_rate_m3_s1"]

            axis.plot(
                np.ma.masked_less_equal(total, 0.0), altitude,
                color="0.12", lw=1.0, ls="--", label="Total", zorder=1,
            )
            for charge, (label, color) in enumerate(
                zip(CHARGE_LABELS, COLORS)
            ):
                axis.plot(
                    np.ma.masked_less_equal(by_charge[:, charge], 0.0),
                    altitude, color=color, lw=1.25, label=label, zorder=2,
                )
            axis.set_xscale("log")
            axis.set_ylim(100, 400)
            axis.grid(True, which="major", color="0.90", lw=0.5)
            axis.set_xlabel(
                rf"{target_labels[row]} ionization rate "
                r"(m$^{-3}$ s$^{-1}$)"
            )
            if row == 0:
                axis.set_title(source_titles[col], fontsize=8)
            axis.text(
                0.025, 0.965, panel_labels[2 * row + col],
                transform=axis.transAxes, ha="left", va="top",
                fontsize=8, fontweight="bold",
            )
            axis.text(
                0.97, 0.965, target_labels[row],
                transform=axis.transAxes, ha="right", va="top",
                fontsize=8, fontweight="bold",
            )

    axes[0, 0].set_ylabel("Altitude (km)")
    axes[1, 0].set_ylabel("Altitude (km)")
    axes[0, 0].legend(
        loc="upper left", bbox_to_anchor=(0.04, 0.90), fontsize=6.5
    )
    fig.suptitle(
        "Target ionization by precipitating hydrogen\n"
        r"100,000 macro particles, 400 km/s, $T=10$ eV, "
        r"$n=5$ cm$^{-3}$",
        fontsize=9,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=600, bbox_inches="tight")
    plt.close(fig)
    print(f"output={args.output.resolve()}")


if __name__ == "__main__":
    main()

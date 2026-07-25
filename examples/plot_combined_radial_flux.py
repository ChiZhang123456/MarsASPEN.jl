"""Plot H ENA-source and H+-source radial flux in a 2 by 3 grid.

Rows identify the injected source state. Columns show downward magnitude,
upward magnitude, and signed outward radial flux from 100 to 400 km.
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

from marsaspen_analysis import load_history_mat, radial_flux_from_mat  # noqa: E402

mpl.rcParams.update({
    "font.family": "Arial",
    "mathtext.fontset": "dejavusans",
    "font.size": 7,
    "axes.linewidth": 0.8,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "legend.frameon": False,
})

SPECIES = (("H ENA", "#4C72B0"), (r"H$^+$", "#C44E52"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("h_ena_mat", type=Path)
    parser.add_argument("hplus_mat", type=Path)
    parser.add_argument(
        "--output",
        type=Path,
        default=REPO / "examples" / "figures" /
        "h_ena_hplus_radial_flux_6panel.png",
    )
    args = parser.parse_args()

    results = [
        radial_flux_from_mat(load_history_mat(args.h_ena_mat)),
        radial_flux_from_mat(load_history_mat(args.hplus_mat)),
    ]
    fig, axes = plt.subplots(
        2, 3, figsize=(7.2, 5.4), sharey=True, constrained_layout=True
    )
    column_titles = (
        "Downward radial flux",
        "Upward radial flux",
        "Net outward radial flux",
    )
    row_labels = ("H ENA source", r"H$^+$ source")
    panel_labels = ("a", "b", "c", "d", "e", "f")

    for row, result in enumerate(results):
        altitude = result["altitude_km"]
        directional = (result["downward"], result["upward"])
        for col in range(3):
            axis = axes[row, col]
            if col < 2:
                for charge, (label, color) in enumerate(SPECIES):
                    axis.plot(
                        np.ma.masked_less_equal(
                            directional[col][:, charge], 0.0
                        ),
                        altitude, color=color, lw=1.25, label=label,
                    )
                axis.set_xscale("log")
                axis.set_xlim(1.0e6, 3.0e12)
                axis.set_xlabel(r"Flux (m$^{-2}$ s$^{-1}$)")
            else:
                for charge, (label, color) in enumerate(SPECIES):
                    axis.plot(
                        result["signed_outward"][:, charge], altitude,
                        color=color, lw=1.25, label=label,
                    )
                axis.axvline(0.0, color="0.35", lw=0.7)
                axis.set_xlim(-2.2e12, 2.2e12)
                axis.ticklabel_format(
                    axis="x", style="sci", scilimits=(0, 0),
                    useMathText=True,
                )
                axis.set_xlabel(
                    r"Signed $F_r$ (m$^{-2}$ s$^{-1}$)"
                )

            axis.set_ylim(100, 400)
            axis.grid(True, which="major", color="0.90", lw=0.5)
            if row == 0:
                axis.set_title(column_titles[col], fontsize=8)
            axis.text(
                0.025, 0.965, panel_labels[3 * row + col],
                transform=axis.transAxes, ha="left", va="top",
                fontsize=8, fontweight="bold",
            )
            if col == 2:
                axis.text(
                    0.97, 0.965, row_labels[row],
                    transform=axis.transAxes, ha="right", va="top",
                    fontsize=7.5, fontweight="bold",
                )

    axes[0, 0].set_ylabel("Altitude (km)")
    axes[1, 0].set_ylabel("Altitude (km)")
    axes[0, 0].legend(loc="best", fontsize=6.5)
    fig.suptitle(
        "Local radial flux of precipitating hydrogen\n"
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

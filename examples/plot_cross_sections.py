"""Plot the collision cross sections used by MarsASPEN.

Rows distinguish the projectile charge state, neutral H ENA or H+. Columns
distinguish the atmospheric target, CO2, O, or N2. Each panel shows the four
channels used by transport: state change, target ionization, Ly-alpha
production, and elastic scattering. Source tables are converted from cm^2 to
m^2 before plotting, matching the SI representation inside MarsASPEN.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

REPO = Path(__file__).resolve().parents[1]
DATA = REPO / "data" / "cross_sections"
OUTPUT = REPO / "examples" / "figures" / "collision_cross_sections.png"

TARGETS = ("CO2", "O", "N2")
CHANNEL_COLORS = ("#C44E52", "#DD8452", "#55A868", "#6B6B6B")
CHANNEL_LABELS = ("State change", "Ionization", r"Ly-$\alpha$", "Elastic")

mpl.rcParams.update({
    "font.family": "Arial",
    "mathtext.fontset": "dejavusans",
    "font.size": 7,
    "axes.linewidth": 0.8,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "legend.frameon": False,
})


def read_table(path: Path) -> np.ndarray:
    """Read a table after its nine comment and header lines."""
    return np.loadtxt(path, comments="#", skiprows=10)


def main() -> None:
    fig, axes = plt.subplots(
        2, 3, figsize=(7.2, 4.7), sharex=True, sharey=True,
        constrained_layout=True,
    )
    projectile_rows = (("H", "H ENA"), ("Hplus", r"H$^+$"))

    for row, (file_prefix, projectile_label) in enumerate(projectile_rows):
        for col, target in enumerate(TARGETS):
            axis = axes[row, col]
            table = read_table(DATA / f"{file_prefix}_{target}_cross_sections.txt")
            energy = table[:, 0]

            # Reorder the source columns into the common transport sequence.
            # H tables: elastic, stripping, Ly-alpha, ionization, Balmer-alpha.
            # H+ tables: charge exchange, ionization, Ly-alpha, elastic,
            # Balmer-alpha. Balmer-alpha is not a transport channel here.
            columns = (2, 4, 3, 1) if file_prefix == "H" else (1, 2, 3, 4)
            for column, label, color in zip(
                columns, CHANNEL_LABELS, CHANNEL_COLORS
            ):
                sigma_m2 = table[:, column] * 1.0e-4
                axis.plot(energy, sigma_m2, color=color, lw=1.1, label=label)

            axis.set_xscale("log")
            axis.set_yscale("log")
            axis.set_xlim(1, 1.0e6)
            axis.set_ylim(1.0e-25, 1.0e-17)
            axis.grid(True, which="major", color="0.90", lw=0.5)
            axis.set_title(f"{projectile_label} + {target}")
            axis.text(
                0.03, 0.95, chr(ord("a") + row * 3 + col),
                transform=axis.transAxes, ha="left", va="top",
                fontweight="bold", fontsize=8,
            )

    for axis in axes[-1, :]:
        axis.set_xlabel("Energy (eV)")
    for axis in axes[:, 0]:
        axis.set_ylabel(r"Cross section (m$^2$)")
    axes[0, 0].legend(ncol=2, fontsize=6, loc="lower left")
    fig.suptitle("MarsASPEN collision cross sections", fontsize=9)
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT, dpi=600, bbox_inches="tight")
    plt.close(fig)
    print(f"output={OUTPUT.resolve()}")


if __name__ == "__main__":
    main()

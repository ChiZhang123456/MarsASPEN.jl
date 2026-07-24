from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

mpl.rcParams["font.family"] = "Arial"
mpl.rcParams["mathtext.fontset"] = "dejavusans"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("csv_file", type=Path)
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args()

    data = np.genfromtxt(args.csv_file, delimiter=",", names=True)
    altitude = data["altitude_center_km"]
    output = args.output or args.csv_file.with_suffix(".png")

    fig, axes = plt.subplots(1, 2, figsize=(12, 7), sharey=True, constrained_layout=True)
    styles = (
        ("state_change", "State change", "#d62728"),
        ("ionization", "Ionization", "#ff7f0e"),
        ("lya", "Ly-alpha", "#2ca02c"),
    )
    for field, label, color in styles:
        axes[0].plot(data[field], altitude, color=color, lw=2, label=label)
    axes[0].set_xlabel("Chemical reaction count per 10 km bin")
    axes[0].set_ylabel("Altitude (km)")
    axes[0].set_xscale("log")
    axes[0].legend()

    axes[1].plot(data["elastic"], altitude, color="#4d4d4d", lw=2, label="Elastic")
    axes[1].plot(
        data["chemical_total"], altitude, color="#1f77b4", lw=2,
        label="All chemical reactions",
    )
    axes[1].set_xlabel("Reaction count per 10 km bin")
    axes[1].set_xscale("log")
    axes[1].legend()
    for ax in axes:
        ax.grid(True, which="both", color="0.85", lw=0.7)
        ax.set_ylim(100, 1000)
    fig.suptitle("MarsASPEN reaction counts, 1,000,000 particles")
    fig.savefig(output, dpi=220)
    plt.close(fig)
    print(f"output={output.resolve()}")


if __name__ == "__main__":
    main()

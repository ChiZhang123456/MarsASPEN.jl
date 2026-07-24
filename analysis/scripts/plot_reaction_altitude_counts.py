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
    parser.add_argument(
        "--reaction",
        choices=("all", "state_change", "ionization", "lya", "elastic"),
        default="all",
    )
    parser.add_argument(
        "--altitude-range", type=float, nargs=2, metavar=("MIN_KM", "MAX_KM")
    )
    args = parser.parse_args()

    data = np.genfromtxt(args.csv_file, delimiter=",", names=True)
    altitude = data["altitude_center_km"]
    output = args.output or args.csv_file.with_suffix(".png")
    bin_width = float(np.median(data["altitude_high_km"] - data["altitude_low_km"]))

    if args.reaction != "all":
        labels = {
            "state_change": "State change",
            "ionization": "Ionization",
            "lya": "Ly-alpha",
            "elastic": "Elastic",
        }
        colors = {
            "state_change": "#d62728",
            "ionization": "#ff7f0e",
            "lya": "#2ca02c",
            "elastic": "#4d4d4d",
        }
        field = args.reaction
        if args.altitude_range is not None:
            fig, ax = plt.subplots(figsize=(7, 7), constrained_layout=True)
            ax.plot(data[field], altitude, color=colors[field], lw=2)
            ax.fill_betweenx(
                altitude, 0, data[field], color=colors[field], alpha=0.18
            )
            ax.set_xlabel(f"{labels[field]} count per {bin_width:g} km bin")
            ax.set_ylabel("Altitude (km)")
            ax.set_ylim(*args.altitude_range)
            ax.grid(True, color="0.85", lw=0.7)
            ax.set_title(
                f"MarsASPEN {labels[field].lower()} counts, "
                "1,000,000 particles"
            )
            fig.savefig(output, dpi=220)
            plt.close(fig)
            print(f"output={output.resolve()}")
            return
        fig, axes = plt.subplots(
            1, 2, figsize=(11, 7), sharey=False, constrained_layout=True
        )
        for ax in axes:
            ax.plot(data[field], altitude, color=colors[field], lw=2)
            ax.fill_betweenx(altitude, 0, data[field], color=colors[field], alpha=0.18)
            ax.set_xlabel(f"{labels[field]} count per {bin_width:g} km bin")
            ax.set_ylabel("Altitude (km)")
            ax.grid(True, color="0.85", lw=0.7)
        axes[0].set_ylim(100, 1000)
        axes[0].set_title("Full altitude range")
        axes[1].set_ylim(100, 250)
        axes[1].set_title("Lower atmosphere detail")
        fig.suptitle(
            f"MarsASPEN {labels[field].lower()} counts, 1,000,000 particles"
        )
        fig.savefig(output, dpi=220)
        plt.close(fig)
        print(f"output={output.resolve()}")
        return

    fig, axes = plt.subplots(1, 2, figsize=(12, 7), sharey=True, constrained_layout=True)
    styles = (
        ("state_change", "State change", "#d62728"),
        ("ionization", "Ionization", "#ff7f0e"),
        ("lya", "Ly-alpha", "#2ca02c"),
    )
    for field, label, color in styles:
        axes[0].plot(data[field], altitude, color=color, lw=2, label=label)
    axes[0].set_xlabel(f"Chemical reaction count per {bin_width:g} km bin")
    axes[0].set_ylabel("Altitude (km)")
    axes[0].set_xscale("log")
    axes[0].legend()

    axes[1].plot(data["elastic"], altitude, color="#4d4d4d", lw=2, label="Elastic")
    axes[1].plot(
        data["chemical_total"], altitude, color="#1f77b4", lw=2,
        label="All chemical reactions",
    )
    axes[1].set_xlabel(f"Reaction count per {bin_width:g} km bin")
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

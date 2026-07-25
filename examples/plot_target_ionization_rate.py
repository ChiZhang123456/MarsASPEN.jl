"""Plot target ionization-rate profiles for H ENA and H+ source simulations."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "analysis"))

from marsaspen_analysis import (  # noqa: E402
    ionization_rate_from_mat,
    load_history_mat,
    mat_string,
)

mpl.rcParams.update({
    "font.family": "Arial",
    "mathtext.fontset": "dejavusans",
    "font.size": 8,
    "axes.linewidth": 0.8,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "legend.frameon": False,
})

COLORS = ("#4C72B0", "#C44E52")
LABELS = ("H ENA contribution", r"H$^+$ contribution")


def source_label(data: dict[str, np.ndarray]) -> str:
    """Return a readable source label stored in a Julia MAT file."""
    value = mat_string(data["initial_species"]).lower()
    return r"H$^+$ source" if "plus" in value else "H ENA source"


def target_labels(data: dict[str, np.ndarray]) -> tuple[str, str]:
    """Return plain-text and Matplotlib labels for the saved target."""
    target = mat_string(data["target_name"]).upper()
    if target == "CO2":
        return "co2", r"CO$_2$"
    if target == "N2":
        return "n2", r"N$_2$"
    return "oxygen", "O"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("h_ena_mat", type=Path)
    parser.add_argument("hplus_mat", type=Path)
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
    )
    args = parser.parse_args()

    datasets = [
        load_history_mat(args.h_ena_mat),
        load_history_mat(args.hplus_mat),
    ]
    results = [ionization_rate_from_mat(data) for data in datasets]
    target_keys = [target_labels(data)[0] for data in datasets]
    if target_keys[0] != target_keys[1]:
        raise ValueError("both MAT files must contain the same target species")
    target_key, target_label = target_labels(datasets[0])

    fig, axes = plt.subplots(
        1, 2, figsize=(5.4, 4.2), sharex=True, sharey=True,
        constrained_layout=True,
    )
    for panel, (axis, data, result) in enumerate(
        zip(axes, datasets, results)
    ):
        altitude = result["altitude_km"]
        by_charge = result["rate_by_charge_m3_s1"]
        total = result["total_rate_m3_s1"]
        axis.plot(
            np.ma.masked_less_equal(total, 0.0), altitude,
            color="0.12", lw=1.1, ls="--", label="Total", zorder=1,
        )
        for charge_index, (label, color) in enumerate(zip(LABELS, COLORS)):
            profile = np.ma.masked_less_equal(by_charge[:, charge_index], 0.0)
            axis.plot(
                profile, altitude, color=color, lw=1.35, label=label, zorder=2
            )
        axis.set_xscale("log")
        axis.set_ylim(80, 600)
        axis.set_xlabel(
            rf"{target_label} ionization rate (m$^{{-3}}$ s$^{{-1}}$)"
        )
        axis.set_title(source_label(data))
        axis.grid(True, which="major", color="0.90", lw=0.55)
        axis.text(
            0.03, 0.97, chr(ord("a") + panel),
            transform=axis.transAxes, ha="left", va="top",
            fontsize=9, fontweight="bold",
        )

    axes[0].set_ylabel("Altitude (km)")
    axes[0].legend(loc="best", fontsize=7)
    fig.suptitle(
        f"{target_label} ionization by precipitating hydrogen\n"
        r"100,000 macro particles, 400 km/s, $T=10$ eV, "
        r"$n=5$ cm$^{-3}$",
        fontsize=9,
    )
    output = args.output or (
        REPO / "examples" / "figures" /
        f"{target_key}_ionization_rate_profiles.png"
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=600, bbox_inches="tight")
    plt.close(fig)
    print(f"output={output.resolve()}")
    for data, result in zip(datasets, results):
        index = int(np.nanargmax(result["total_rate_m3_s1"]))
        print(
            f"{source_label(data)}: peak="
            f"{result['total_rate_m3_s1'][index]:.9g} m^-3 s^-1 "
            f"at {result['altitude_km'][index]:.1f} km"
        )


if __name__ == "__main__":
    main()

"""Plot Ly-alpha VER, radiative energy, and optically thin limb brightness."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "analysis"))

from marsaspen_analysis import load_history_mat, lya_profiles_from_mat  # noqa: E402

mpl.rcParams.update({
    "font.family": "Arial",
    "mathtext.fontset": "dejavusans",
    "font.size": 7,
    "axes.linewidth": 0.8,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "legend.frameon": False,
})

SPECIES = (("H ENA contribution", "#4C72B0"),
           (r"H$^+$ contribution", "#C44E52"))
XLIMITS = (
    (1.0e2, 1.0e9),
    (1.0e-16, 1.0e-8),
    (1.0e-2, 1.0e5),
)


def plot_source_profile(result, source_label: str, output: Path) -> None:
    """Save one source simulation as a three-panel altitude profile figure."""
    altitude = result["altitude_km"]
    profiles = (
        result["volume_emission_rate_by_charge_photons_m3_s1"],
        result["radiative_energy_rate_by_charge_w_m3"],
        result["limb_brightness_by_charge_rayleigh"],
    )
    totals = (
        result["total_volume_emission_rate_photons_m3_s1"],
        result["total_radiative_energy_rate_w_m3"],
        result["total_limb_brightness_rayleigh"],
    )
    titles = (
        "Ly-alpha volume emission rate",
        "Ly-alpha radiative energy",
        "Optically thin limb brightness",
    )
    xlabels = (
        r"VER (photons m$^{-3}$ s$^{-1}$)",
        r"Radiative energy rate (W m$^{-3}$)",
        "Brightness (R)",
    )
    fig, axes = plt.subplots(
        1, 3, figsize=(7.2, 2.8), sharey=True, constrained_layout=True
    )
    for col, axis in enumerate(axes):
        axis.plot(
            np.ma.masked_less_equal(totals[col], 0.0), altitude,
            color="0.12", lw=1.0, ls="--", label="Total", zorder=3,
        )
        for charge, (label, color) in enumerate(SPECIES):
            axis.plot(
                np.ma.masked_less_equal(profiles[col][:, charge], 0.0),
                altitude, color=color, lw=1.25, label=label, zorder=2,
            )
        axis.set_xscale("log")
        axis.set_xlim(*XLIMITS[col])
        axis.set_ylim(100, 400)
        axis.set_xlabel(xlabels[col])
        axis.set_title(titles[col], fontsize=8)
        axis.grid(True, which="major", color="0.90", lw=0.5)
        axis.text(
            0.025, 0.965, "abc"[col], transform=axis.transAxes,
            ha="left", va="top", fontsize=8, fontweight="bold",
        )
    axes[0].set_ylabel("Altitude (km)")
    axes[0].legend(loc="upper left", bbox_to_anchor=(0.04, 0.90), fontsize=6.5)
    fig.suptitle(
        f"H Ly-alpha profiles, {source_label}\n"
        r"100,000 macro particles, 400 km/s, $T=10$ eV, "
        r"$n=5$ cm$^{-3}$",
        fontsize=9,
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=600, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("h_ena_mat", type=Path)
    parser.add_argument("hplus_mat", type=Path)
    parser.add_argument(
        "--output",
        type=Path,
        default=REPO / "examples" / "figures" /
        "h_ena_hplus_lya_profiles_6panel.png",
    )
    args = parser.parse_args()

    results = [
        lya_profiles_from_mat(load_history_mat(args.h_ena_mat)),
        lya_profiles_from_mat(load_history_mat(args.hplus_mat)),
    ]
    fig, axes = plt.subplots(
        2, 3, figsize=(7.2, 5.4), sharex="col", sharey=True,
        constrained_layout=True,
    )
    column_titles = (
        "Ly-alpha volume emission rate",
        "Ly-alpha radiative energy",
        "Optically thin limb brightness",
    )
    row_labels = ("H ENA source", r"H$^+$ source")
    panel_labels = ("a", "b", "c", "d", "e", "f")

    for row, result in enumerate(results):
        altitude = result["altitude_km"]
        profiles = (
            result["volume_emission_rate_by_charge_photons_m3_s1"],
            result["radiative_energy_rate_by_charge_w_m3"],
            result["limb_brightness_by_charge_rayleigh"],
        )
        totals = (
            result["total_volume_emission_rate_photons_m3_s1"],
            result["total_radiative_energy_rate_w_m3"],
            result["total_limb_brightness_rayleigh"],
        )
        xlabels = (
            r"VER (photons m$^{-3}$ s$^{-1}$)",
            r"Radiative energy rate (W m$^{-3}$)",
            "Brightness (R)",
        )

        for col in range(3):
            axis = axes[row, col]
            axis.plot(
                np.ma.masked_less_equal(totals[col], 0.0), altitude,
                color="0.12", lw=1.0, ls="--", label="Total", zorder=3,
            )
            for charge, (label, color) in enumerate(SPECIES):
                axis.plot(
                    np.ma.masked_less_equal(profiles[col][:, charge], 0.0),
                    altitude, color=color, lw=1.25, label=label, zorder=2,
                )
            axis.set_xscale("log")
            axis.set_xlim(*XLIMITS[col])
            axis.set_ylim(100, 400)
            if row == 1:
                axis.set_xlabel(xlabels[col])
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
    axes[0, 0].legend(
        loc="upper left", bbox_to_anchor=(0.04, 0.90), fontsize=6.5
    )
    fig.suptitle(
        "H Ly-alpha production by precipitating hydrogen\n"
        r"100,000 macro particles, 400 km/s, $T=10$ eV, "
        r"$n=5$ cm$^{-3}$",
        fontsize=9,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=600, bbox_inches="tight")
    plt.close(fig)
    plot_source_profile(
        results[0], "H ENA source",
        args.output.with_name("h_ena_lya_profiles_3panel.png"),
    )
    plot_source_profile(
        results[1], r"H$^+$ source",
        args.output.with_name("hplus_lya_profiles_3panel.png"),
    )
    print(f"output={args.output.resolve()}")
    for label, result in zip(row_labels, results):
        peak = int(np.nanargmax(
            result["total_volume_emission_rate_photons_m3_s1"]
        ))
        print(
            f"{label}: peak_VER="
            f"{result['total_volume_emission_rate_photons_m3_s1'][peak]:.9g} "
            f"photons m^-3 s^-1 at {result['altitude_km'][peak]:.1f} km; "
            f"max_limb={np.nanmax(result['total_limb_brightness_rayleigh']):.9g} R"
        )


if __name__ == "__main__":
    main()

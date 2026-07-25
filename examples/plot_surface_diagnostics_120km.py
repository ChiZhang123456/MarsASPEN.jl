"""Plot 120 km point-source ionization, Ly-alpha, and energy footprints."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm

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
})


def positive_limits(arrays: list[np.ndarray]) -> tuple[float, float]:
    """Return four-decade shared limits that resolve the footprint core."""
    positive = np.concatenate([a[np.isfinite(a) & (a > 0)] for a in arrays])
    if positive.size == 0:
        return 1.0, 10.0
    high = 10.0 ** np.ceil(np.log10(np.max(positive)))
    return high / 1.0e4, high


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("h_ena_mat", type=Path)
    parser.add_argument("hplus_mat", type=Path)
    parser.add_argument(
        "--output", type=Path,
        default=REPO / "examples" / "figures" /
        "point_source_surface_diagnostics_120km_6panel.png",
    )
    args = parser.parse_args()
    data = [
        load_history_mat(args.h_ena_mat),
        load_history_mat(args.hplus_mat),
    ]
    keys = (
        "total_ionization_rate_m3_s1",
        "lya_volume_emission_rate_photons_m3_s1",
        "total_energy_transfer_w_m3",
    )
    titles = (
        "O + CO$_2$ ionization rate",
        "H Ly-alpha volume emission rate",
        "Projectile energy transfer rate",
    )
    colorbar_labels = (
        r"Ionization rate (m$^{-3}$ s$^{-1}$)",
        r"VER (photons m$^{-3}$ s$^{-1}$)",
        r"Energy transfer (W m$^{-3}$)",
    )
    limits = [
        positive_limits([np.asarray(d[key], dtype=float) for d in data])
        for key in keys
    ]

    fig, axes = plt.subplots(
        2, 3, figsize=(7.2, 4.7), sharex=True, sharey=True,
        constrained_layout=True,
    )
    row_labels = ("H ENA source", r"H$^+$ source")
    panel_labels = "abcdef"
    for row, source in enumerate(data):
        lon_edges = np.ravel(source["longitude_edges_deg"])
        lat_edges = np.ravel(source["latitude_edges_deg"])
        for col, key in enumerate(keys):
            axis = axes[row, col]
            values = np.asarray(source[key], dtype=float)
            masked = np.ma.masked_less_equal(values.T, 0.0)
            mesh = axis.pcolormesh(
                lon_edges, lat_edges, masked, shading="flat",
                cmap="turbo", norm=LogNorm(*limits[col]),
                rasterized=True,
            )
            axis.set_aspect("equal")
            axis.set_xlim(lon_edges[0], lon_edges[-1])
            axis.set_ylim(lat_edges[0], lat_edges[-1])
            axis.set_xlabel("MSO longitude (deg)")
            if col == 0:
                axis.set_ylabel("MSO latitude (deg)")
            if row == 0:
                axis.set_title(titles[col], fontsize=8)
            axis.text(
                0.025, 0.975, panel_labels[3 * row + col],
                transform=axis.transAxes, ha="left", va="top",
                fontsize=8, fontweight="bold",
            )
            if col == 2:
                axis.text(
                    0.975, 0.975, row_labels[row],
                    transform=axis.transAxes, ha="right", va="top",
                    fontsize=7.5, fontweight="bold",
                )
            colorbar = fig.colorbar(mesh, ax=axis, pad=0.02, fraction=0.046)
            colorbar.set_label(colorbar_labels[col])

    fig.suptitle(
        "Point-source precipitation footprint at 120 km\n"
        r"100,000 macro particles, 400 km/s, $T=10$ eV, "
        r"$n=5$ cm$^{-3}$",
        fontsize=9,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=600, bbox_inches="tight")
    plt.close(fig)
    print(f"output={args.output.resolve()}")
    for key, (low, high) in zip(keys, limits):
        print(f"{key}: clim=({low:.6g}, {high:.6g})")


if __name__ == "__main__":
    main()

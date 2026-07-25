"""Plot the weighted H+ energy distribution at one selected altitude.

This diagnostic reads the compact altitude-energy MAT output produced by
``run_h_ena_100000_monte_carlo.jl``. The plotted quantity is density-weighted
track length divided by the energy-bin width, evaluated in the one-kilometer
altitude bin containing the requested altitude.
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

from marsaspen_analysis.io import load_history_mat  # noqa: E402

mpl.rcParams["font.family"] = "Arial"
mpl.rcParams["mathtext.fontset"] = "dejavusans"


def vector(data: dict[str, np.ndarray], key: str) -> np.ndarray:
    """Return one MAT vector as a one-dimensional floating-point array."""
    return np.asarray(data[key], dtype=float).squeeze()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("mat_file", type=Path)
    parser.add_argument("--altitude", type=float, default=550.0)
    parser.add_argument(
        "--rebin",
        type=int,
        default=5,
        help="Number of native energy bins combined for the smooth curve.",
    )
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args()

    data = load_history_mat(args.mat_file)
    altitude_edges = vector(data, "altitude_edges_km")
    energy_edges = vector(data, "energy_edges_ev")
    track_length = np.asarray(data["weighted_track_length_m4"], dtype=float)

    expected_shape = (
        altitude_edges.size - 1,
        energy_edges.size - 1,
        2,
    )
    if track_length.shape == expected_shape[::-1]:
        track_length = track_length.transpose(2, 1, 0)
    if track_length.shape != expected_shape:
        raise ValueError(
            f"Unexpected histogram shape {track_length.shape}; "
            f"expected {expected_shape}."
        )

    altitude_index = np.searchsorted(
        altitude_edges, args.altitude, side="right"
    ) - 1
    if not 0 <= altitude_index < altitude_edges.size - 1:
        raise ValueError("Requested altitude is outside the MAT altitude grid.")

    energy_centers = 0.5 * (energy_edges[:-1] + energy_edges[1:])
    energy_widths = np.diff(energy_edges)
    # Charge index 1 is H+. Division by dE gives a differential distribution.
    distribution = track_length[altitude_index, :, 1] / energy_widths
    positive = distribution > 0
    if not np.any(positive):
        raise ValueError("No H+ contribution exists in the selected altitude bin.")

    if args.rebin < 1:
        raise ValueError("--rebin must be at least one.")
    usable_bins = distribution.size // args.rebin * args.rebin
    rebinned_distribution = (
        track_length[altitude_index, :usable_bins, 1]
        .reshape(-1, args.rebin)
        .sum(axis=1)
        / (
            energy_edges[args.rebin:usable_bins + 1:args.rebin]
            - energy_edges[:usable_bins:args.rebin]
        )
    )
    rebinned_energy = 0.5 * (
        energy_edges[:usable_bins:args.rebin]
        + energy_edges[args.rebin:usable_bins + 1:args.rebin]
    )
    peak_index = int(np.argmax(rebinned_distribution))
    peak_energy = rebinned_energy[peak_index]
    altitude_low = altitude_edges[altitude_index]
    altitude_high = altitude_edges[altitude_index + 1]

    fig, axis = plt.subplots(figsize=(8.5, 5.5), constrained_layout=True)
    axis.step(
        energy_centers,
        distribution,
        where="mid",
        color="0.65",
        linewidth=0.8,
        label=f"Native bins, {energy_widths[0]:g} eV",
    )
    axis.step(
        rebinned_energy,
        rebinned_distribution,
        where="mid",
        color="#d62728",
        linewidth=1.8,
        label=(
            f"Rebinned, "
            f"{args.rebin * energy_widths[0]:g} eV"
        ),
    )
    axis.axvline(
        peak_energy,
        color="black",
        linestyle="--",
        linewidth=1.0,
        label=f"Rebinned peak = {peak_energy:.1f} eV",
    )
    axis.set_xlim(10, 1200)
    axis.set_ylim(bottom=0)
    axis.set_xlabel("H+ energy (eV)")
    axis.set_ylabel(
        r"Density-weighted track length / $\Delta E$ "
        r"(m$^{-2}$ eV$^{-1}$)"
    )
    axis.set_title(
        f"H+ energy distribution at {altitude_low:.0f}–"
        f"{altitude_high:.0f} km\n"
        "100,000 initially neutral H ENA particles"
    )
    axis.grid(True, color="0.85", linewidth=0.7)
    axis.legend()

    output = args.output or args.mat_file.with_name(
        f"hplus_energy_{altitude_low:.0f}_{altitude_high:.0f}km.png"
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=220)
    plt.close(fig)

    print(f"altitude_bin_km={altitude_low:.1f},{altitude_high:.1f}")
    print(f"nonzero_energy_bins={np.count_nonzero(positive)}")
    print(f"rebinned_energy_width_ev={args.rebin * energy_widths[0]:.1f}")
    print(f"rebinned_peak_energy_ev={peak_energy:.1f}")
    print(f"weighted_track_length_m2={track_length[altitude_index, :, 1].sum():.9g}")
    print(f"output={output.resolve()}")


if __name__ == "__main__":
    main()

"""Plot the weighted H+ energy distribution at one selected altitude.

This diagnostic reads the compact altitude-energy MAT output produced by
``run_h_ena_100000_monte_carlo.jl``. The plotted quantity is the direct sum of
particle density weights, evaluated in the one-kilometer altitude bin
containing the requested altitude.

Native logarithmic bins and a rebinned curve are both shown. Rebinning changes
only this diagnostic display. It does not rerun transport or modify MAT data.
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

    # Read grid edges and the altitude × energy × charge histogram.
    data = load_history_mat(args.mat_file)
    altitude_edges = vector(data, "altitude_edges_km")
    energy_edges = vector(data, "energy_edges_ev")
    density_weight_sum = np.asarray(data["density_weight_sum_m3"], dtype=float)

    expected_shape = (
        altitude_edges.size - 1,
        energy_edges.size - 1,
        2,
    )
    # h5py exposes Julia array dimensions in reverse order for MAT v7.3.
    if density_weight_sum.shape == expected_shape[::-1]:
        density_weight_sum = density_weight_sum.transpose(2, 1, 0)
    if density_weight_sum.shape != expected_shape:
        raise ValueError(
            f"Unexpected histogram shape {density_weight_sum.shape}; "
            f"expected {expected_shape}."
        )

    # Locate the one-kilometer layer containing the requested altitude.
    altitude_index = np.searchsorted(
        altitude_edges, args.altitude, side="right"
    ) - 1
    if not 0 <= altitude_index < altitude_edges.size - 1:
        raise ValueError("Requested altitude is outside the MAT altitude grid.")

    energy_centers = 0.5 * (energy_edges[:-1] + energy_edges[1:])
    # Charge index 1 is H+. Each bin is already the requested direct sum.
    distribution = density_weight_sum[altitude_index, :, 1]
    positive = distribution > 0
    if not np.any(positive):
        raise ValueError("No H+ contribution exists in the selected altitude bin.")

    # Sum neighboring bins because the saved value is integrated per energy
    # bin, not divided by the bin width.
    if args.rebin < 1:
        raise ValueError("--rebin must be at least one.")
    usable_bins = distribution.size // args.rebin * args.rebin
    rebinned_distribution = (
        density_weight_sum[altitude_index, :usable_bins, 1]
        .reshape(-1, args.rebin)
        .sum(axis=1)
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
        label="Native logarithmic bins",
    )
    axis.step(
        rebinned_energy,
        rebinned_distribution,
        where="mid",
        color="#d62728",
        linewidth=1.8,
        label=f"Combined every {args.rebin} bins",
    )
    axis.axvline(
        peak_energy,
        color="black",
        linestyle="--",
        linewidth=1.0,
        label=f"Rebinned peak = {peak_energy:.1f} eV",
    )
    axis.set_xscale("log")
    axis.set_xlim(1, 10_000)
    axis.set_ylim(bottom=0)
    axis.set_xlabel("H+ energy (eV)")
    axis.set_ylabel(r"Sum of particle density weights (m$^{-3}$)")
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
    print(f"combined_native_bins={args.rebin}")
    print(f"rebinned_peak_energy_ev={peak_energy:.1f}")
    print(
        "density_weight_sum_m3="
        f"{density_weight_sum[altitude_index, :, 1].sum():.9g}"
    )
    print(f"output={output.resolve()}")


if __name__ == "__main__":
    main()

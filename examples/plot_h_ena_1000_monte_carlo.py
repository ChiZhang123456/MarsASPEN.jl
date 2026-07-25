"""Plot the weighted altitude-energy distribution from the 1,000 H ENA run.

The Julia example stores density-weighted track length in each altitude,
energy, and charge-state bin. This script plots H ENA and H+ separately. H+
appears when the initially neutral projectile undergoes electron stripping.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "analysis"))

from marsaspen_analysis.io import load_history_mat  # noqa: E402

# Apply the project figure convention to every non-mathematical text object.
mpl.rcParams["font.family"] = "Arial"
mpl.rcParams["mathtext.fontset"] = "dejavusans"


def one_dimensional(data: dict[str, np.ndarray], key: str) -> np.ndarray:
    """Read a MATLAB row or column vector as a NumPy one-dimensional array."""
    return np.asarray(data[key], dtype=float).squeeze()


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Plot the weighted 1,000-particle H ENA example."
    )
    parser.add_argument("mat_file", type=Path)
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args()

    # The shared reader supports both ordinary MAT files and MAT v7.3/HDF5.
    data = load_history_mat(args.mat_file)
    altitude_edges = one_dimensional(data, "altitude_edges_km")
    energy_edges = one_dimensional(data, "energy_edges_ev")
    weighted_track_length = np.asarray(
        data["weighted_track_length_m4"], dtype=float
    )
    expected_shape = (
        altitude_edges.size - 1,
        energy_edges.size - 1,
        2,
    )
    if weighted_track_length.shape == expected_shape[::-1]:
        # HDF5 exposes Julia/MATLAB dimensions in reverse order.
        weighted_track_length = weighted_track_length.transpose(2, 1, 0)
    if weighted_track_length.shape != expected_shape:
        raise ValueError(
            "Unexpected weighted_track_length_m4 shape: "
            f"{weighted_track_length.shape}, expected {expected_shape}."
        )

    # A shared logarithmic color scale makes the H ENA and H+ panels directly
    # comparable. Empty bins are masked rather than assigned an artificial
    # floor.
    positive = weighted_track_length[weighted_track_length > 0]
    if positive.size == 0:
        raise ValueError("The MAT file contains no positive weighted bins.")
    color_norm = LogNorm(vmin=np.percentile(positive, 2), vmax=positive.max())

    fig, axes = plt.subplots(
        1, 2, figsize=(12, 6), sharex=True, sharey=True,
        constrained_layout=True,
    )
    species_names = ("H ENA", r"H$^+$")
    image = None
    for charge_index, (axis, species) in enumerate(zip(axes, species_names)):
        values = np.ma.masked_less_equal(
            weighted_track_length[:, :, charge_index], 0
        )
        image = axis.pcolormesh(
            energy_edges,
            altitude_edges,
            values,
            shading="auto",
            cmap="viridis",
            norm=color_norm,
        )
        axis.set_title(species)
        axis.set_xlabel("Energy (eV)")
        axis.set_xlim(10, 1000)
        axis.set_ylim(80, 300)
        axis.grid(False)

    axes[0].set_ylabel("Altitude (km)")
    colorbar = fig.colorbar(image, ax=axes, pad=0.02)
    colorbar.set_label(r"Density-weighted track length (m$^{-2}$)")
    fig.suptitle(
        "1,000 H ENA particles, 600 km, 400 km/s, "
        r"$T=10$ eV, $n=1$ cm$^{-3}$"
    )

    output = args.output or args.mat_file.with_suffix(".png")
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=220)
    plt.close(fig)
    print(f"output={output.resolve()}")


if __name__ == "__main__":
    main()

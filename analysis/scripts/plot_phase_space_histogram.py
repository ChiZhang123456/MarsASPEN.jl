"""Plot path-length-weighted H ENA and H+ altitude-energy histograms."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import h5py
from matplotlib.colors import LogNorm

mpl.rcParams["font.family"] = "Arial"
mpl.rcParams["mathtext.fontset"] = "dejavusans"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("mat_file", type=Path)
    parser.add_argument("--output", type=Path, default=None)
    parser.add_argument("--altitude-range", type=float, nargs=2, default=(80, 600))
    args = parser.parse_args()

    with h5py.File(args.mat_file, "r") as data:
        altitude_edges = np.asarray(data["altitude_edges_km"], dtype=float).ravel()
        energy_edges = np.asarray(data["energy_edges_ev"], dtype=float).ravel()
        # Julia MAT v7.3 arrays appear in reversed dimension order through h5py.
        path_length_km = (
            np.asarray(data["path_length_m"], dtype=float).transpose(2, 1, 0)
            / 1000.0
        )
        n_particles = int(np.asarray(data["n_particles"]).item())
        particle_weight = float(np.asarray(data["particle_weight"]).item())

    positive = path_length_km[path_length_km > 0]
    norm = LogNorm(
        vmin=max(np.percentile(positive, 1), positive.min()),
        vmax=positive.max(),
    )
    fig, axes = plt.subplots(
        1, 2, figsize=(12, 7), sharex=True, sharey=True, constrained_layout=True
    )
    titles = ("H ENA", r"H$^+$")
    for charge, (ax, title) in enumerate(zip(axes, titles)):
        mesh = ax.pcolormesh(
            energy_edges,
            altitude_edges,
            path_length_km[:, :, charge],
            shading="flat",
            cmap="viridis",
            norm=norm,
        )
        ax.set_xscale("log")
        ax.set_xlabel("Energy (eV)")
        ax.set_title(title)
        ax.grid(False)
    axes[0].set_ylabel("Altitude (km)")
    axes[0].set_ylim(*args.altitude_range)
    cbar = fig.colorbar(mesh, ax=axes, pad=0.02)
    cbar.set_label(r"Weighted trajectory path length, $\sum w_i ds$ (km)")
    fig.suptitle(
        f"MarsASPEN altitude-energy distribution, {n_particles:,} particles, "
        f"$w_i$ = {particle_weight:g}"
    )
    output = args.output or args.mat_file.with_suffix(".png")
    fig.savefig(output, dpi=220)
    plt.close(fig)
    print(f"output={output.resolve()}")


if __name__ == "__main__":
    main()

"""Plot upward and downward H ENA and H+ differential number flux."""

from __future__ import annotations

import argparse
from pathlib import Path

import h5py
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm

mpl.rcParams["font.family"] = "Arial"
mpl.rcParams["mathtext.fontset"] = "dejavusans"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("mat_file", type=Path)
    parser.add_argument("--output", type=Path, default=None)
    parser.add_argument("--altitude-range", type=float, nargs=2, default=(100, 300))
    parser.add_argument("--energy-range", type=float, nargs=2, default=(10, 1400))
    args = parser.parse_args()

    with h5py.File(args.mat_file, "r") as data:
        altitude = np.asarray(data["altitude_surfaces_km"], dtype=float).ravel()
        energy_edges = np.asarray(data["energy_edges_ev"], dtype=float).ravel()
        # Julia dimensions are reversed in the HDF5 representation.
        flux = np.asarray(data["flux_m2_s"], dtype=float).transpose(3, 2, 1, 0)
        n_particles = int(np.asarray(data["n_particles"]).item())
        density_cm3 = (
            float(np.asarray(data["initial_number_density_m3"]).item()) / 1e6
        )
        temperature_ev = float(np.asarray(data["initial_temperature_ev"]).item())
        bulk_velocity = np.asarray(data["initial_bulk_velocity_m_s"], dtype=float).ravel()

    energy_center = 0.5 * (energy_edges[:-1] + energy_edges[1:])
    energy_width = np.diff(energy_edges)
    differential_flux = flux / energy_width[None, :, None, None]

    selected = (
        (altitude >= args.altitude_range[0])
        & (altitude <= args.altitude_range[1])
    )
    plot_values = differential_flux[selected]
    positive = plot_values[plot_values > 0]
    norm = LogNorm(
        vmin=max(np.percentile(positive, 0.5), positive.min()),
        vmax=positive.max(),
    )

    fig, axes = plt.subplots(
        2, 2, figsize=(12, 10), sharex=True, sharey=True,
        constrained_layout=True,
    )
    panel_specs = (
        (1, 0, r"Downward H$^+$"),
        (1, 1, r"Upward H$^+$"),
        (0, 0, "Downward H ENA"),
        (0, 1, "Upward H ENA"),
    )
    mesh = None
    for ax, (charge, direction, title) in zip(axes.ravel(), panel_specs):
        mesh = ax.pcolormesh(
            energy_center,
            altitude,
            differential_flux[:, :, charge, direction],
            shading="nearest",
            cmap="viridis",
            norm=norm,
        )
        ax.set_title(title)
        ax.set_xlim(*args.energy_range)
        ax.set_ylim(*args.altitude_range)
    for ax in axes[-1, :]:
        ax.set_xlabel("Energy (eV)")
    for ax in axes[:, 0]:
        ax.set_ylabel("Altitude (km)")
    cbar = fig.colorbar(mesh, ax=axes, pad=0.02)
    cbar.set_label(
        r"Differential number flux (m$^{-2}$ s$^{-1}$ eV$^{-1}$)"
    )
    fig.suptitle(
        "MarsASPEN solar-wind proton transport, "
        f"{n_particles:,} particles, n = {density_cm3:g} cm$^{{-3}}$, "
        f"kT = {temperature_ev:g} eV, "
        f"V = [{bulk_velocity[0]/1000:g}, 0, 0] km/s"
    )
    output = args.output or args.mat_file.with_suffix(".png")
    fig.savefig(output, dpi=220)
    plt.close(fig)
    print(f"output={output.resolve()}")


if __name__ == "__main__":
    main()

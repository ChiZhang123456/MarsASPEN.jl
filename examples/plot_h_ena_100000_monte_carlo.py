"""Plot a 100,000-particle density-weight altitude-energy distribution.

Each color bin is the direct sum of particle density weights for crossings in
that altitude and energy bin. No path length or energy-width normalization is
applied. The reader accepts either the initially neutral H ENA example or the
initially ionized H+ example and constructs the title and color limits from the
metadata stored by Julia.

Important interpretation: the current histogram combines upward and downward
crossings. A particle that returns through the same altitude is counted again.
Consequently, the plotted quantity is useful for comparing the simulated
populations, but it is not a unique-particle count or a local density estimate.
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


def matlab_string(data: dict[str, np.ndarray], key: str) -> str:
    """Decode a string written by MAT.jl in either MAT or MAT v7.3 format."""
    value = np.asarray(data[key]).squeeze()
    if value.dtype.kind in "US":
        return "".join(value.reshape(-1).astype(str))
    if value.dtype.kind in "ui":
        return "".join(chr(int(code)) for code in value.reshape(-1) if code)
    return str(value)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Plot a weighted MarsASPEN altitude-energy MAT output."
    )
    parser.add_argument("mat_file", type=Path)
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args()

    # The shared reader supports both ordinary MAT files and MAT v7.3/HDF5.
    data = load_history_mat(args.mat_file)
    altitude_edges = one_dimensional(data, "altitude_edges_km")
    energy_edges = one_dimensional(data, "energy_edges_ev")
    density_weight_sum = np.asarray(
        data["density_weight_sum_m3"], dtype=float
    )
    expected_shape = (
        altitude_edges.size - 1,
        energy_edges.size - 1,
        2,
    )
    if density_weight_sum.shape == expected_shape[::-1]:
        # HDF5 exposes Julia/MATLAB dimensions in reverse order.
        density_weight_sum = density_weight_sum.transpose(2, 1, 0)
    if density_weight_sum.shape != expected_shape:
        raise ValueError(
            "Unexpected density_weight_sum_m3 shape: "
            f"{density_weight_sum.shape}, expected {expected_shape}."
        )

    # Positive values are required by LogNorm. Empty bins are masked below
    # rather than replaced with a small artificial number.
    if not np.any(density_weight_sum > 0):
        raise ValueError("The MAT file contains no positive weighted bins.")

    # Read the physical source metadata rather than hard-coding a title. For
    # the H ENA source, the requested display ranges remain 1e1 to 1e6 m^-3
    # for H ENA and 1e1 to 1e5 m^-3 for the smaller H+ product. For the 5 cm^-3
    # proton source, H+ is the dominant population, so its upper limit is
    # raised to 1e7 m^-3 while H ENA retains the 1e6 m^-3 limit.
    initial_species = matlab_string(data, "initial_species")
    n_particles = int(one_dimensional(data, "n_particles"))
    initial_altitude = float(one_dimensional(data, "initial_altitude_km"))
    temperature = float(one_dimensional(data, "physical_temperature_ev"))
    source_density = float(one_dimensional(data, "source_number_density_m3"))
    bulk_speed = np.linalg.norm(
        one_dimensional(data, "initial_bulk_velocity_m_s")
    ) / 1000.0
    initially_ionized = "plus" in initial_species.lower()
    color_norms = (
        LogNorm(vmin=1.0e1, vmax=1.0e6, clip=True),
        LogNorm(vmin=1.0e1, vmax=1.0e7 if initially_ionized else 1.0e5, clip=True),
    )

    fig, axes = plt.subplots(
        1, 2, figsize=(12, 6), sharex=True, sharey=True,
        constrained_layout=True,
    )
    species_names = ("H ENA", r"H$^+$")
    for charge_index, (axis, species, color_norm) in enumerate(
        zip(axes, species_names, color_norms)
    ):
        # Charge index zero is neutral H ENA and index one is H+. Masking zeros
        # keeps genuinely empty regions white on the logarithmic color scale.
        values = np.ma.masked_less_equal(
            density_weight_sum[:, :, charge_index], 0
        )
        image = axis.pcolormesh(
            energy_edges,
            altitude_edges,
            values,
            shading="auto",
            cmap="turbo",
            norm=color_norm,
        )
        axis.set_title(species)
        axis.set_xlabel("Energy (eV)")
        axis.set_xscale("log")
        axis.set_xlim(10, 3_000)
        axis.set_ylim(100, 600)
        axis.grid(False)
        colorbar = fig.colorbar(image, ax=axis, pad=0.02)
        colorbar.set_label(r"Sum of particle density weights (m$^{-3}$)")

    axes[0].set_ylabel("Altitude (km)")
    source_label = r"H$^+$" if initially_ionized else "H ENA"
    fig.suptitle(
        f"{n_particles:,} {source_label} particles, "
        f"{initial_altitude:g} km, {bulk_speed:g} km/s, "
        rf"$T={temperature:g}$ eV, $n={source_density / 1.0e6:g}$ cm$^{{-3}}$"
    )

    output = args.output or args.mat_file.with_suffix(".png")
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=220)
    plt.close(fig)
    print(f"output={output.resolve()}")


if __name__ == "__main__":
    main()

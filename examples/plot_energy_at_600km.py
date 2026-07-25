"""Diagnose the weighted particle energy distribution just below 600 km."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import ncx2

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "analysis"))

from marsaspen_analysis.io import load_history_mat  # noqa: E402

mpl.rcParams["font.family"] = "Arial"
mpl.rcParams["mathtext.fontset"] = "dejavusans"


def vector(data: dict[str, np.ndarray], key: str) -> np.ndarray:
    """Return a MAT vector as a one-dimensional floating-point array."""
    return np.asarray(data[key], dtype=float).squeeze()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("mat_file", type=Path)
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args()

    data = load_history_mat(args.mat_file)
    altitude_edges = vector(data, "altitude_edges_km")
    energy_edges = vector(data, "energy_edges_ev")
    histogram = np.asarray(data["density_weight_sum_m3"], dtype=float)
    expected_shape = (altitude_edges.size - 1, energy_edges.size - 1, 2)
    if histogram.shape == expected_shape[::-1]:
        histogram = histogram.transpose(2, 1, 0)
    if histogram.shape != expected_shape:
        raise ValueError(
            f"Unexpected histogram shape {histogram.shape}, expected {expected_shape}."
        )

    # The saved altitude values are bin edges. The highest crossing surface is
    # 599.5 km, represented by the final 599 to 600 km row.
    altitude_index = altitude_edges.size - 2
    h_ena = histogram[altitude_index, :, 0]
    hplus = histogram[altitude_index, :, 1]

    temperature_ev = float(vector(data, "physical_temperature_ev"))
    source_density = float(vector(data, "source_number_density_m3"))
    speed = np.linalg.norm(vector(data, "initial_bulk_velocity_m_s"))
    proton_mass_kg = 1.67262192369e-27
    elementary_charge_c = 1.602176634e-19
    bulk_energy_ev = 0.5 * proton_mass_kg * speed**2 / elementary_charge_c

    # For a three-dimensional drifting Maxwellian, 2E/T follows a noncentral
    # chi-square distribution with three degrees of freedom and
    # noncentrality 2 E_bulk/T. This is the injection distribution expected
    # for downward particles before atmospheric processing.
    scaled_edges = 2.0 * energy_edges / temperature_ev
    noncentrality = 2.0 * bulk_energy_ev / temperature_ev
    expected_injection = source_density * np.diff(
        ncx2.cdf(scaled_edges, df=3, nc=noncentrality)
    )

    centers = np.sqrt(energy_edges[:-1] * energy_edges[1:])
    fig, axis = plt.subplots(figsize=(9, 5.8), constrained_layout=True)
    axis.step(
        centers, h_ena, where="mid", linewidth=1.8, color="#0072B2",
        label="Simulated H ENA, upward + downward crossings",
    )
    axis.step(
        centers, hplus, where="mid", linewidth=1.5, color="#D55E00",
        label=r"Simulated H$^+$, upward + downward crossings",
    )
    axis.step(
        centers, expected_injection, where="mid", linewidth=1.5,
        color="black", linestyle="--",
        label="Expected downward injection, 10 eV drifting Maxwellian",
    )
    axis.axvline(
        bulk_energy_ev, color="0.35", linestyle=":", linewidth=1.2,
        label=f"400 km/s bulk energy = {bulk_energy_ev:.1f} eV",
    )
    axis.set_xscale("log")
    axis.set_yscale("log")
    axis.set_xlim(1, 10_000)
    axis.set_ylim(bottom=1.0e-6)
    axis.set_xlabel("Energy (eV)")
    axis.set_ylabel(r"Sum of particle density weights (m$^{-3}$)")
    axis.set_title(
        "Energy distribution at the 599.5 km crossing surface\n"
        "100,000 particles, native 30 logarithmic energy bins"
    )
    axis.grid(True, which="both", color="0.88", linewidth=0.6)
    axis.legend(fontsize=9)

    output = args.output or args.mat_file.with_name("energy_distribution_600km.png")
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=220)
    plt.close(fig)

    print(f"altitude_surface_km=599.5")
    print(f"bulk_energy_ev={bulk_energy_ev:.9g}")
    print(f"h_ena_density_weight_sum_m3={h_ena.sum():.9g}")
    print(f"hplus_density_weight_sum_m3={hplus.sum():.9g}")
    print(f"expected_injection_sum_m3={expected_injection.sum():.9g}")
    print(f"output={output.resolve()}")


if __name__ == "__main__":
    main()

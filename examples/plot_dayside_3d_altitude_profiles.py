"""Plot global spherical-mean profiles from the full 3D Monte Carlo output."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "analysis"))
ALTITUDE_MIN_KM = 100.0
ALTITUDE_MAX_KM = 300.0

from marsaspen_analysis import load_history_mat  # noqa: E402

mpl.rcParams.update({
    "font.family": "Arial",
    "mathtext.fontset": "dejavusans",
    "font.size": 7,
    "axes.linewidth": 0.8,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "legend.frameon": False,
})


def spherical_area_mean(
    values: np.ndarray,
    latitude_edges_deg: np.ndarray,
    longitude_centers_deg: np.ndarray,
) -> np.ndarray:
    """Return a spherical-area-weighted horizontal mean at every altitude."""
    latitude_weight = np.diff(
        np.sin(np.deg2rad(latitude_edges_deg))
    )
    longitude_weight = np.ones(longitude_centers_deg.size, dtype=float)
    weights = latitude_weight[:, None] * longitude_weight[None, :]
    return np.sum(values * weights[None, :, :], axis=(1, 2)) / np.sum(weights)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "prefix", nargs="?", type=Path,
        default=REPO / "examples" / "output" /
        "dayside_hplus_10000000_3d",
    )
    parser.add_argument(
        "--output", type=Path,
        default=REPO / "examples" / "figures" /
        "dayside_hplus_10000000_3d_altitude_profiles_6panel.png",
    )
    args = parser.parse_args()
    prefix = str(args.prefix)
    moments = load_history_mat(Path(prefix + "_moments.mat"))
    reactions = load_history_mat(Path(prefix + "_reactions.mat"))
    energy = load_history_mat(Path(prefix + "_energy.mat"))

    altitude = np.ravel(moments["altitude_centers_km"])
    latitude_edges = np.ravel(moments["latitude_edges_deg"])
    longitude_centers = np.ravel(moments["longitude_centers_deg"])
    density = np.asarray(moments["number_density_by_charge_m3"], dtype=float)
    downward = np.asarray(
        moments["downward_radial_flux_by_charge_m2_s"], dtype=float
    )
    upward = np.asarray(
        moments["upward_radial_flux_by_charge_m2_s"], dtype=float
    )
    ionization = np.asarray(
        reactions["ionization_rate_by_target_m3_s1"], dtype=float
    )
    lya = np.asarray(
        reactions["total_lya_volume_emission_rate_photons_m3_s1"],
        dtype=float,
    )
    deposition = np.asarray(energy["total_energy_transfer_w_m3"], dtype=float)

    mean = lambda values: spherical_area_mean(  # noqa: E731
        values, latitude_edges, longitude_centers
    )
    profiles = {
        "hena_density": mean(density[0]),
        "hplus_density": mean(density[1]),
        "hplus_down": mean(downward[1]),
        "hplus_up": mean(upward[1]),
        "hena_down": mean(downward[0]),
        "hena_up": mean(upward[0]),
        "co2_ionization": mean(ionization[0]),
        "o_ionization": mean(ionization[1]),
        "lya": mean(lya),
        "deposition": mean(deposition),
    }

    panels = (
        (
            "H$^+$ and H-ENA number density",
            r"Number density (m$^{-3}$)",
            (
                (profiles["hplus_density"], "H$^+$", "#C44E52", "-"),
                (profiles["hena_density"], "H-ENA", "#4C72B0", "-"),
            ),
        ),
        (
            "H$^+$ radial flux",
            r"Radial flux (m$^{-2}$ s$^{-1}$)",
            (
                (profiles["hplus_down"], "Downward", "#C44E52", "-"),
                (profiles["hplus_up"], "Upward", "#C44E52", "--"),
            ),
        ),
        (
            "H-ENA radial flux",
            r"Radial flux (m$^{-2}$ s$^{-1}$)",
            (
                (profiles["hena_down"], "Downward", "#4C72B0", "-"),
                (profiles["hena_up"], "Upward", "#4C72B0", "--"),
            ),
        ),
        (
            "Target ionization rate",
            r"Ionization rate (m$^{-3}$ s$^{-1}$)",
            (
                (profiles["o_ionization"], "O", "#2CA02C", "-"),
                (profiles["co2_ionization"], "CO$_2$", "#9467BD", "-"),
            ),
        ),
        (
            "H Ly-alpha volume emission rate",
            r"VER (photons m$^{-3}$ s$^{-1}$)",
            ((profiles["lya"], "H Ly-alpha", "#D18F00", "-"),),
        ),
        (
            "Energy deposition rate",
            r"Energy deposition rate (W m$^{-3}$)",
            ((profiles["deposition"], "Total", "0.15", "-"),),
        ),
    )

    fig, axes = plt.subplots(
        2, 3, figsize=(7.2, 5.0), sharey=True, constrained_layout=True
    )
    for panel, (axis, (title, xlabel, curves)) in enumerate(
        zip(axes.flat, panels)
    ):
        axis.set_xscale("log")
        for values, label, color, linestyle in curves:
            axis.plot(
                np.ma.masked_less_equal(values, 0.0), altitude,
                color=color, lw=1.2, ls=linestyle, label=label,
            )
        axis.set(
            ylim=(ALTITUDE_MIN_KM, ALTITUDE_MAX_KM),
            title=title, xlabel=xlabel,
        )
        if panel % 3 == 0:
            axis.set_ylabel("Altitude (km)")
        axis.grid(True, which="major", color="0.90", lw=0.5)
        axis.legend(loc="upper left", fontsize=6.5)
        axis.text(
            0.025, 0.975, "abcdef"[panel],
            transform=axis.transAxes, ha="left", va="top",
            fontsize=8, fontweight="bold",
        )
    fig.suptitle(
        "Global spherical-mean profiles from uniform-dayside H$^+$ injection\n"
        f"{ALTITUDE_MIN_KM:.0f} to {ALTITUDE_MAX_KM:.0f} km, "
        r"10,000,000 particles, $\mathbf{U}=(-400,0,0)$ km/s, "
        r"$T=10$ eV, $n=5$ cm$^{-3}$, "
        r"5$^\circ\times$5$^\circ\times1$ km grid",
        fontsize=9,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=600, bbox_inches="tight")
    plt.close(fig)
    print(f"output={args.output.resolve()}")
    displayed = (
        (altitude >= ALTITUDE_MIN_KM) &
        (altitude <= ALTITUDE_MAX_KM)
    )
    for name, profile in profiles.items():
        if np.any(np.isfinite(profile[displayed])):
            displayed_indices = np.flatnonzero(displayed)
            index = displayed_indices[
                int(np.nanargmax(np.abs(profile[displayed])))
            ]
            print(
                f"{name}: global_peak={profile[index]:.9g} "
                f"at {altitude[index]:.1f} km"
            )


if __name__ == "__main__":
    main()

"""Plot area-weighted altitude profiles from the full 3D Monte Carlo output."""

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
    dayside_only: bool,
) -> np.ndarray:
    """Return a spherical-area-weighted horizontal mean at every altitude."""
    latitude_weight = np.diff(
        np.sin(np.deg2rad(latitude_edges_deg))
    )
    longitude_mso = (
        (longitude_centers_deg + 180.0) % 360.0
    ) - 180.0
    longitude_mask = (
        np.abs(longitude_mso) < 90.0
        if dayside_only else np.ones(longitude_mso.size, dtype=bool)
    )
    weights = latitude_weight[:, None] * longitude_mask[None, :]
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
        "dayside_hplus_10000000_3d_altitude_profiles_8panel.png",
    )
    args = parser.parse_args()
    prefix = str(args.prefix)
    moments = load_history_mat(Path(prefix + "_moments.mat"))
    reactions = load_history_mat(Path(prefix + "_reactions.mat"))
    energy = load_history_mat(Path(prefix + "_energy.mat"))

    altitude = np.ravel(moments["altitude_centers_km"])
    latitude_edges = np.ravel(moments["latitude_edges_deg"])
    longitude_centers = np.ravel(moments["longitude_centers_deg"])
    volumes = [
        np.asarray(moments["total_number_density_m3"], dtype=float),
        np.asarray(moments["total_flux_m2_s"], dtype=float),
        np.asarray(moments["downward_radial_flux_m2_s"], dtype=float),
        np.asarray(moments["upward_radial_flux_m2_s"], dtype=float),
        np.asarray(moments["signed_radial_flux_m2_s"], dtype=float),
        np.asarray(reactions["reaction_rate_by_channel_m3_s1"], dtype=float)[1],
        np.asarray(
            reactions["total_lya_volume_emission_rate_photons_m3_s1"],
            dtype=float,
        ),
        np.asarray(energy["total_energy_transfer_w_m3"], dtype=float),
    ]
    global_profiles = [
        spherical_area_mean(
            values, latitude_edges, longitude_centers, False
        )
        for values in volumes
    ]
    dayside_profiles = [
        spherical_area_mean(
            values, latitude_edges, longitude_centers, True
        )
        for values in volumes
    ]
    titles = (
        "H + H$^+$ number density",
        "Total scalar flux",
        "Downward radial flux",
        "Upward radial flux",
        "Signed outward radial flux",
        "Total target ionization rate",
        "H Ly-alpha volume emission rate",
        "Projectile energy transfer rate",
    )
    xlabels = (
        r"Density (m$^{-3}$)",
        r"Flux (m$^{-2}$ s$^{-1}$)",
        r"Downward flux (m$^{-2}$ s$^{-1}$)",
        r"Upward flux (m$^{-2}$ s$^{-1}$)",
        r"Signed flux (m$^{-2}$ s$^{-1}$)",
        r"Ionization rate (m$^{-3}$ s$^{-1}$)",
        r"VER (photons m$^{-3}$ s$^{-1}$)",
        r"Energy transfer (W m$^{-3}$)",
    )

    fig, axes = plt.subplots(
        2, 4, figsize=(10.0, 5.3), sharey=True, constrained_layout=True
    )
    for panel, (axis, title, xlabel, global_profile, dayside_profile) in enumerate(
        zip(
            axes.flat, titles, xlabels,
            global_profiles, dayside_profiles,
        )
    ):
        if panel == 4:
            axis.set_xscale("symlog", linthresh=1.0e7)
            global_values = global_profile
            dayside_values = dayside_profile
        else:
            axis.set_xscale("log")
            global_values = np.ma.masked_less_equal(global_profile, 0.0)
            dayside_values = np.ma.masked_less_equal(dayside_profile, 0.0)
        axis.plot(
            global_values, altitude, color="0.15", lw=1.1,
            label="Global spherical mean",
        )
        axis.plot(
            dayside_values, altitude, color="#C44E52", lw=1.2,
            label="Dayside mean",
        )
        axis.set(
            ylim=(ALTITUDE_MIN_KM, ALTITUDE_MAX_KM),
            title=title, xlabel=xlabel,
        )
        if panel % 4 == 0:
            axis.set_ylabel("Altitude (km)")
        axis.grid(True, which="major", color="0.90", lw=0.5)
        axis.text(
            0.025, 0.975, "abcdefgh"[panel],
            transform=axis.transAxes, ha="left", va="top",
            fontsize=8, fontweight="bold",
        )
    axes[0, 0].legend(loc="upper left", bbox_to_anchor=(0.08, 0.92), fontsize=6.5)
    fig.suptitle(
        "Uniform-dayside H$^+$ Monte Carlo altitude profiles\n"
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
    for title, profile in zip(titles, dayside_profiles):
        if np.any(np.isfinite(profile[displayed])):
            displayed_indices = np.flatnonzero(displayed)
            index = displayed_indices[
                int(np.nanargmax(np.abs(profile[displayed])))
            ]
            print(
                f"{title}: dayside_peak={profile[index]:.9g} "
                f"at {altitude[index]:.1f} km"
            )


if __name__ == "__main__":
    main()
